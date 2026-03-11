#' Compute clustering performance metrics
#'
#' This function evaluates the quality of the clustering results by comparing the 
#' estimated community assignments against the ground truth derived from the 
#' adjacency matrix using Normalized Mutual Information (NMI), Adjusted Mutual 
#' Information (AMI), and Adjusted Rand Index (ARI).
#'
#' @param clus_member A numeric or integer vector of estimated cluster assignments.
#' @param A A square adjacency matrix (ground truth). 
#'
#' @return A named list containing three numeric scores:
#' \itemize{
#'   \item \code{NMI}: Normalized Mutual Information
#'   \item \code{AMI}: Adjusted Mutual Information
#'   \item \code{ARI}: Adjusted Rand Index
#' }
#' 
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom aricode NMI AMI ARI
#' @export
#'
ComputePerformance <- function(clus_member, A){
  # first, find the clusters in the adjacency matrix
  graph <- graph_from_adjacency_matrix(A, mode="undirected")
  clusters <- components(graph)$membership
  
  NMI <- NMI(clus_member,clusters)
  AMI <- AMI(clus_member,clusters)
  ARI <- ARI(clus_member,clusters)
  
  score <- list(NMI,AMI,ARI)
  names(score) <- c("NMI","AMI","ARI")
  return(score)
}
