# code l1 spectral clustering
library(glmnet)

l1_spectralclustering <- function(A,k){
  m <- ncol(A)
  
  # 1st step: SVD
  svd <- eigen(A) 
  U <- svd$vectors
  
  # 2nd step: solve the lasso problem to find the first community
  w <- U[1,(k+1):m]
  W <- t(U[-1,(k+1):m])
  
  lambda <- cv.glmnet(W,w)
  lambda <- lambda$lambda.min
  v <- glmnet(W,w,lambda=lambda)
  v <- c(1,v)
  
  # 3rd step: find the other communities using deflation
  
}