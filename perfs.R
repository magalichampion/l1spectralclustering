# code for testing the performances of the algorithms
score <- Compute_NMIscore(comm=results$comm,A=data$A,method="l1")
score 

Compute_NMIscore <- function(comm,A,method){
  # comm: community indices
  # A: adjacency matrix of the true graph
  # method: l1-spetral, spectral,...
  real <- RealComm(A)
  clus_real <- ClustVec(real)
  clus_real <- cbind(c(1:length(clus_real)),clus_real)
  if (method=="l1"){
    clus_est <- ClustVec(comm)
  } else if (method=="spectral"){
    clus_est <- comm
  }
  clus_est <- cbind(c(1:length(clus_est)),clus_est)
  score <- NMI(clus_est,clus_real)
  score <- score$value
  return(score)
}

RealComm <- function(A){
  mat <- A+diag(dim(A)[2])
  cluster <- c()
  for (i in (1:length(unique(colSums(mat))))){
    s <- unique(colSums(mat))[i]
    clus <- length(which(colSums(mat)==s))/s
    cluster <- c(cluster,rep(s,clus))
  }
  perfect <- matrix(0,sum(cluster),length(cluster))
  for(i in 1:length(cluster)){
    perfect[,i]<- mat[,cumsum(cluster)[i]]
  }
  return(perfect)
}

ClustVec <- function(comm){
  if (!is.null(ncol(comm))){
    new <- comm%*%c(1:dim(comm)[2])
  } else {
    new <- comm
  }
  return(new)
}
