#' Create a simulated perturbed graph and the associated adjacency matrix
#'
#' @description 
#' Generates a synthetic graph with a ground-truth block-diagonal structure. 
#' The function allows for controlled perturbations by removing edges from 
#' within clusters (p_inside) and adding noisy edges between clusters (p_outside).
#'
#' @param k Integer. The number of ground-truth clusters.
#' @param n Integer. The total number of nodes in the graph.
#' @param p List. A list containing two numeric values:
#'   \itemize{
#'     \item \code{p_inside}: Probability of removing an existing edge within a cluster.
#'     \item \code{p_outside}: Probability of adding a noise edge between clusters.
#'   }
#' @param print.plot Logical. If \code{TRUE}, displays side-by-side plots of the ideal vs. perturbed graphs.
#' @param ClustersLength Numeric vector. Optional manual specification of cluster sizes. 
#' If \code{NULL}, sizes are generated randomly such that they sum to \code{n}.
#' @param rngSEED Integer. To use as seed.
#' 
#' @return A list containing:
#' \item{A}{The ideal (ground-truth) adjacency matrix.}
#' \item{A_hat}{The perturbed/noisy adjacency matrix used for clustering benchmarks.}
#' \item{ClustersLength}{The size of each generated cluster.}
#'
#' @importFrom Matrix bdiag
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph
#' 
#' @export
#'
#' @examples
#' perturb <- list(p_inside = 0.1, p_outside = 0.05)
#' data <- CreateDataSet(k = 3, n = 30, p = perturb)
#' 
CreateDataSet <- function(k, n, p, print.plot=TRUE, ClustersLength = NULL, rngSEED){
  p1 <- p$p_inside
  p2 <- p$p_outside
  
  if (k>n){
    stop("The number of clusters is larger than the number of nodes. Please choose a smaller number of clusters.")
  } 
  
  MaxLength <- n - 2*(k-1) 
  MinLength <- 2
  if (MinLength!=MaxLength){
    Length <- sample(MinLength:MaxLength)
  } else {
    Length <- 2
  }
  
  if (is.null(ClustersLength)){
    ClustersLength <- c()
    Length_tmp <- Length
    for (i in (1:(k-1))){
      Length_tmp <- Length_tmp[Length_tmp<(n-sum(ClustersLength)+2*(i-k)+1)]
      if (length(Length_tmp)>0){
        clustersLength <- Length_tmp[1]
        ClustersLength <- c(ClustersLength,clustersLength)
        Length_tmp <- Length_tmp[-1]
      } else {
        ClustersLength <- c(ClustersLength,2)
      }
    }
    ClustersLength <- c(ClustersLength,(n-sum(ClustersLength)))
    ClustersLength <- sort(ClustersLength)
  }
  
  print(paste0(c("There are",k,"clusters of size ",ClustersLength,"."), collapse=" "))
  
  A <- matrix(1,ncol=ClustersLength[1],nrow=ClustersLength[1])
  for (i in 2:k){
    A <- bdiag(A,matrix(1,ClustersLength[i],ClustersLength[i]));
  }
  A <- matrix(A,ncol(A),nrow(A))
  
  A <- A-diag(nrow(A)) # adjacency matrix
  
  # perturbed versions of the graph  
  A_perturbed <- matrix(0,n,n)
  edges_inside = edges_outside = 0
  bern_inside = bern_outside <- c()
  for (i in (1:k)){
    bern_inside <- rbinom(ClustersLength[i]*(ClustersLength[i]-1)/2,size=1,prob=1-p1)
    if (i>1){
      bern_outside <- rbinom(sum(ClustersLength[1:(i-1)])*ClustersLength[i],size=1,prob=p2)
    }
    
    B <- matrix(0,ClustersLength[i],ClustersLength[i])
    B[upper.tri(B, diag = FALSE)] <- bern_inside
    
    if (i>1){
      A_perturbed[(sum(ClustersLength[1:(i-1)])+1):sum(ClustersLength[1:i]),(sum(ClustersLength[1:(i-1)])+1):sum(ClustersLength[1:i])] <- B
      A_perturbed[(1:sum(ClustersLength[1:(i-1)])),(sum(ClustersLength[1:(i-1)])+1):(sum(ClustersLength[1:i]))] <- bern_outside
    } else {
      A_perturbed[(1:ClustersLength[1]),(1:ClustersLength[1])] <- B
    }
    edges_inside <- edges_inside + length(which(bern_inside==0))
    edges_outside <- edges_outside + length(which(bern_outside==1))
  }
  A_perturbed <- A_perturbed+t(A_perturbed)
  
  print(paste0("On the ",length(which(A==1))/2," existing edges, ",edges_inside," were removed and ",edges_outside," were added."))
  
  par(mfrow=c(1,2))
  graph_hat <- graph_from_adjacency_matrix(A_perturbed,mode="undirected")
  graph <- graph_from_adjacency_matrix(A,mode="undirected")
  plot(graph,main="Real graph")
  plot(graph_hat,main="Perturbed graph")
  
  data <- list(A=A, A_hat = A_perturbed, ClustersLength=ClustersLength)
}
