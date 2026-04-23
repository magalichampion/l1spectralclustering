#' Run cluster analysis for a single data set
#'
#' This is the execution function for the simulation study. It accepts a graph 
#' object, applies a specified clustering algorithm, measures execution time, 
#' and calculates performance metrics.
#'
#' @param graph_obj A list containing at least:
#'   \itemize{
#'     \item \code{A_hat}: The observed noisy adjacency matrix.
#'     \item \code{A}: The ground-truth adjacency matrix (used for performance scoring).
#'     \item \code{ClustersLength}: A vector indicating the size of each ground-truth cluster.
#'   }
#' @param method A character string specifying the algorithm to use. 
#'   Options include:
#'   \itemize{
#'     \item \code{"l1Spectral"}: l1-penalized clustering.
#'     \item \code{"Spectral"}: Standard spectral clustering.
#'     \item \code{"regSpectral"}: Regularized spectral clustering.
#'     \item \code{"robustSpectral"}: Robust spectral clustering.
#'     \item \code{"VGAE}: Variational Graph AutoEncoder.
#'     \item \code{"ST_l1Spectral"}: Self-tuned version of l1-spectral clustering.
#'     \item \code{"Spectrum"}: Self-tuned spectral clustering.
#'     \item \code{"Hybrid"}: Hybrid clustering algorithm.
#'     \item \code{"MCL"}: Markov Cluster Algorithm.
#'   }
#'
#' @return A named list containing:
#'   \item{clusters}{A vector of cluster assignments (or \code{NA} if the method failed).}
#'   \item{NMI}{Normalized Mutual Information score.}
#'   \item{AMI}{Adjusted Mutual Information score.}
#'   \item{ARI}{Adjusted Rand Index score.}
#'   \item{Time}{Execution time in seconds.}
#'
#' @export

RunSimulations <- function(graph_obj, method = c("l1Spectral", "Spectral", "regSpectral","robustSpectral", "VGAE", "ST_l1Spectral", "Spectrum", "Hybrid", "MCL")) {
  # Start timer
  t1 <- Sys.time()
  
  # Run chosen method
  clusters <- switch(method,
                "l1Spectral" = l1Spectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "Spectral" = Spectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "regSpectral" = regSpectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "robustSpectral" = robustSpectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "VGAE" = VGAE(graph_obj$A_hat,length(graph_obj$ClustersLength)),
                "ST_l1Spectral" = l1Spectral(graph_obj$A_hat),
                "Spectrum" = Spectrum(graph_obj$A_hat),
                "Hybrid" = Hybrid(graph_obj$A_hat),
                "MCL" = MCL(graph_obj$A_hat)
  )
  
  # End timer
  t2 <- Sys.time()
  tdiff <- as.numeric(difftime(t2, t1, units = "secs"))
  
  # Check if clusters is a single NA, NULL, or an empty vector
  if (is.null(clusters) || (length(clusters) == 1 && is.na(clusters))) {
    return(list(
      clusters = NA,
      NMI = NA,
      AMI = NA,
      ARI = NA,
      Time = tdiff  
    ))
  } else {
    # Compute the performance
    scores <- ComputePerformance(clusters, graph_obj$A)
  
    # Return structured list
    return(list(
      clusters = clusters,
      NMI = scores$NMI,
      AMI = scores$AMI,
      ARI = scores$ARI,
      Time = tdiff
    ))
  }
}

#' l1-spectral clustering 
#'
#' The l1-spectral clustering algorithm from the l1spectral package
#'
#' @param A A square adjacency matrix.
#' @param k The number of clusters to find (NULL by default).
#'
#' @return A list of clustering results returned by \code{l1spectral::l1_spectralclustering}. 
#'   Typically includes cluster assignments.
#' 
#' @importFrom l1spectral l1_spectralclustering
#' @importFrom purrr possibly
#' @export
#'
l1Spectral <- function(A,k = NULL){
  
  run_l1Spectral <- function(mat, centers) {
    Structure <- l1spectral::FindStructure(A = mat)
    k_hat <- l1spectral::FindNbrClusters(A = mat, structure  = Structure, k = centers, k_max = 25)
    results_l1Spectral <- l1spectral::l1_spectralclustering(A=mat,pen="thresholdedLS",k=k_hat$nbr_clusters_total ,k_max=25)
    if (!is.null(ncol(results_l1Spectral$comm))){
      clusters <- results_l1Spectral$comm%*%c(1:ncol(results_l1Spectral$comm))
    } else {
      clusters <- results_l1Spectral$comm
    }
    clusters <- as.vector(clusters)
    return(clusters)
  }
  
  # If it fails, return NA so RunSimulations can handle it
  safe_l1Spectral <- purrr::possibly(run_l1Spectral, otherwise = NA)
  
  results <- safe_l1Spectral(A,k)
  
  return(results)
}          

#' Spectral clustering algorithm
#'
#' The basic spectral algorithm from the kernlab package.
#'
#' @param A An adjacency matrix.
#' @param k The number of clusters to find.
#' @return A numeric vector of cluster assignments, or NA if the algorithm fails.
#' @importFrom kernlab specc
#' @importFrom purrr possibly

Spectral <- function(A, k) {
  
  run_specc <- function(mat, centers) {
    res <- kernlab::specc(mat, centers = centers)
    return(as.vector(res@.Data))
  }
  
  # Wrap with 'possibly' to handle errors (e.g., singular matrices)
  # We set the 'otherwise' to NA
  safe_specc <- purrr::possibly(run_specc, otherwise = NA)
  
  results <- safe_specc(A, centers = k)
  
  return(results)
}

#' Regularized spectral clustering algorithm
#'
#' The regularized spectral clustering algorithm from the greed package.
#'
#' @param A An adjacency matrix.
#' @param k The number of clusters to find.
#' @return A numeric vector of cluster assignments.
#' @importFrom greed spectral
#' @export
regSpectral <- function(A, k) {
  results <- greed::spectral(A, K = k)
  
  # Ensure the results match the number of nodes
  # Isolated nodes are often dropped, we reassign them to cluster 0 
  if (length(results) < ncol(A)) {
    full_results <- rep(0, ncol(A))
    full_results[1:length(results)] <- results
    results <- full_results
  }
  
  return(results)
}

#' Robust spectral clustering algorithm
#'
#' The robust spectral clustering algorithm (Python code)
#'
#' @param A An adjacency matrix.
#' @param k The number of clusters to find.
#' @return A numeric vector of cluster assignments.
#'
#' @importFrom reticulate source_python
#' @export
robustSpectral <- function(A, k) {
  python_path <- here::here("python", "rsc_code.py")
  reticulate::source_python(python_path)
  
  A_mat <- matrix(as.numeric(as.matrix(A)),nrow=nrow(A))
  clusters <- run_rsc(A_mat, as.integer(k))
  clusters <- clusters + 1 # to rename the clusters starting from 1
  clusters <- as.vector(clusters)
  return(clusters)
}

#' Variational Graph Auto-Encoder (VGAE) 
#'
#' This function implements a Variational Graph Auto-Encoder using a GCN encoder 
#' and an inner-product decoder. It maps an adjacency matrix into a latent space 
#' and performs clustering on the resulting embeddings (Python code)
#'
#' @param A An adjacency matrix.
#' @param k The number of clusters to find.
#' @return A numeric vector of cluster assignments.
#' 
#' @export
#' @importFrom here here
#' @importFrom reticulate py_run_string import py_run_file
VGAE <- function(A,k){
  
  vgae_path <- here("python", "vgae_pytorch")
  py_run_string(sprintf("import sys; sys.path.append('%s')", vgae_path))
  
  # define parameters
  n_nodes <- nrow(A) 
  py_run_string(sprintf("
import args
args.input_dim = %d
args.hidden1_dim = 32
args.hidden2_dim = 16
", as.integer(nrow(A))))
  
  # process the adj matrix
  # D^-0.5 * (A + I) * D^-0.5
  A_tilde <- A + diag(n_nodes)
  D_inv_sqrt <- diag(rowSums(A_tilde)^-0.5)
  A_norm_mat <- D_inv_sqrt %*% A_tilde %*% D_inv_sqrt
  
  # convert to Torch Tensors
  torch <- import("torch")
  adj_norm  <- torch$tensor(as.matrix(A_norm_mat), dtype=torch$float32)
  
  torch_f <- import("torch.nn.functional")
  
  adj_label <- torch$tensor(as.matrix(A + diag(n_nodes)), dtype=torch$float32)
  features  <- torch$eye(as.integer(n_nodes), dtype=torch$float32)
  
  # initialize and train
  model_lib <- import("model")
  vgae_inst <- model_lib$VGAE(adj_norm)
  optimizer <- torch$optim$Adam(vgae_inst$parameters(), lr=0.01)
  
  pos_weight <- (n_nodes^2 - sum(A)) / sum(A)
  weight_vec <- rep(1, n_nodes^2)
  adj_flat <- as.vector(as.matrix(A + diag(n_nodes)))
  edge_indices <- which(adj_flat == 1)
  weight_vec[edge_indices] <- pos_weight
  weight_tensor <- torch$tensor(weight_vec, dtype = torch$float32)
  
  total_entries <- n_nodes * n_nodes
  non_edges <- total_entries - sum(A)
  norm_val <- total_entries / (non_edges * 2)
  
  for (epoch in 1:200) {
    vgae_inst$train()
    optimizer$zero_grad()
    
    # Forward pass
    A_pred <- vgae_inst(features)
    A_pred_flat <- A_pred$reshape(-1L)
    
    # Loss Calculation (BCE + KL)
    adj_label_flat <- adj_label$contiguous()$reshape(-1L)
    
    # 3. Now run the loss
    log_lik <- norm_val * torch_f$binary_cross_entropy(A_pred_flat, adj_label_flat, weight = weight_tensor)
    kl <- 0.5 / n_nodes * torch$sum(1 + 2*vgae_inst$logstd - vgae_inst$mean^2 - torch$exp(vgae_inst$logstd)^2, 1L)$mean()
    
    loss <- log_lik - kl
    loss$backward()
    optimizer$step()
    
    if(epoch %% 20 == 0) cat("Epoch:", epoch, "Loss:", loss$item(), "\n")
  }
  
  vgae_inst$eval()
  embeddings <- vgae_inst$encode(features)$detach()$numpy()
  
  # Use KMeans to find the clusters
  clusters <- kmeans(embeddings, centers = k, nstart = 25)$cluster
}

#' Spectral Clustering using the Spectrum algorithm
#'
#' The spectrum algorithm from the \code{Spectrum} package.
#'
#' @param A An adjacency matrix.
#' @return A numeric vector of cluster assignments. 
#'
#' @importFrom Spectrum Spectrum
#' @export
#'
Spectrum <- function(A){
  results <- Spectrum::Spectrum(A,kerneltype = "stsc")
  clusters <- c(results$assignments)
}

#' Hybrid clustering algorithm (ICL-based)
#'
#' The hybrid clustering algorithm from the \code{greed} package.
#'
#' @param A An adjacency matrix.
#' @return A numeric vector of cluster assignments. 
#' @importFrom greed greed
#' @export
Hybrid <- function(A){
  results <- greed::greed(A)
  clusters <- c(results@cl)
}

#' Markov Clustering algorithm (MCL)
#'
#' The MCL algorithm from the MCL package
#'
#' @param A An adjacency matrix.
#' @return A numeric vector of cluster assignments.
#' @importFrom MCL mcl
#' @importFrom purrr possibly
#' @export
MCL <- function(A){
  run_mcl <- function(mat) {
    mcl_results <- MCL::mcl(mat, addLoops = TRUE)
    return(mcl_results$Cluster)
  }
  
  # If it fails, return NA so RunSimulations can handle it
  safe_mcl <- purrr::possibly(run_mcl, otherwise = NA)
  
  results <- safe_mcl(A)
  return(results)
}
