# code l1 spectral clustering
library(glmnet)

l1_spectralclustering <- function(A,k){
  m <- ncol(A)
  
  # 1 com
  res <- find_com(A,k,m)
  com <- res$v
  
  for (i in (1:k)){
    res <- find_com(res$Adef,k,m) 
    com <- cbind(com,res$v)
  }
}

find_com <- function(A,k,m){
  # 1st step: SVD
  svd <- eigen(A) 
  U <- svd$vectors
  
  # 2nd step: solve the lasso problem to find the first community
  w <- U[1,(k+1):m]
  W <- t(U[-1,(k+1):m])
  
  lassosol <- cv.glmnet(W,w)
  if (lassosol$nzero[which(lassosol$lambda==lassosol$lambda.min)]==0){
    Firstnnzero <- lassosol$lambda[which(lassosol$nzero>0)]
    lambda_opt <- Firstnnzero[1]
  } else {
    lambda_opt <- lassosol$lambda.min
  }
  sol <- glmnet(W,w,lambda=lambda_opt)
  sol <- sol$beta
  
  v <- c(1,as.matrix(sol))
  
  # 3rd step: find the other communities using deflation
  A_def <- A - t(t(v)) %*% t(v)
  
  return(list(A_def,v))
}