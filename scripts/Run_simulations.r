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
#'     \item \code{"ST_l1Spectral"}: Self-tuned version of l1-spectral clustering.
#'     \item \code{"ST_Spectral"}: Self-tuned spectral clustering.
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

RunSimulations <- function(graph_obj, method = c("l1Spectral", "Spectral", "regSpectral","robustSpectral","ST_l1spectral", "ST_Spectral", "Hybrid", "MCL")) {
  # Start timer
  t1 <- Sys.time()
  
  # Run chosen method
  clusters <- switch(method,
                "l1Spectral" = l1Spectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "Spectral" = Spectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "regSpectral" = regSpectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "robustSpectral" = robustSpectral(graph_obj$A_hat, length(graph_obj$ClustersLength)),
                "ST_l1Spectral" = l1Spectral(graph_obj$A_hat),
                "ST_Spectral" = ST_Spectral(graph_obj$A_hat),
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
#' @export
#'
l1Spectral <- function(A,k = NULL){
  run_l1Spectral <- function(mat, centers) {
    results_l1Spectral <- l1spectral::l1_spectralclustering(A=mat,pen="thresholdedLS",k=centers, k_max = 50)
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
  
  # Call the Python function
  clusters <- run_rsc(as.matrix(A), k)
  
  return(clusters)
}

# ST spectral clustering
ST_Spectral <- function(A){
  python_script <- here("python", "stsc.py")
  
  source_python(python_script)
  
  sol <- self_tuning_spectral_clustering(A)
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