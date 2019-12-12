# code l1 spectral clustering
library(glmnet)

l1_spectralclustering <- function(A,k,index){
  # inputs :
  # - A : matrice d'adjacence
  # - k : nombre de communautés
  # - index : vecteur indiquant 1 noeud par communauté
  
  # nombre de noeuds
  m <- ncol(A)
  
  # 1st step: travail sur la 1ère communauté
  res <- find_com(A,k,m,index,iteration=1)
  com <- res$vector
  
  # 2nd step: loop sur les communautés
  for (i in (2:k)){
    res <- find_com(res$A,k,m,index,iteration=i) 
    com <- cbind(com,res$vector)
  }
  
  com[which(com<0.5)] <- 0
  com[which(com>0.5)] <- 1
}

find_com <- function(A,k,m,index,iteration){
  # main code pour trouver les indices de chaque communauté
  
  # 1st step: SVD sur A 
  svd <- eigen(A) 
  eigenvalues <- sort(svd$values,index.return=TRUE)
  eigenvectors <- svd$vectors[,eigenvalues$ix]
  U <- t(eigenvectors[,1:(m-k+iteration-1)])
  
  # 2nd step: solve the lasso problem to find the first community
  w <- U[,index[iteration]]
  W <- U[,-index[iteration]]
  
  #lassosol <- cv.glmnet(W,-w)
  # if (lassosol$nzero[which(lassosol$lambda==lassosol$lambda.min)]==0){
  #   Firstnnzero <- lassosol$lambda[which(lassosol$nzero>0)]
  #   lambda_opt <- Firstnnzero[1]
  # } else {
  #   lambda_opt <- lassosol$lambda.min
  # }
 
  #lambda_opt <- lassosol$lambda[which(lassosol$cvm>min(lassosol$cvup))]
  #cvtop<- min(lassosol$cvm)+5*(max(lassosol$cvm)-min(lassosol$cvm))/100
  #plot(lassosol)
  #abline(h=cvtop)
  #lambda_opt <- lassosol$lambda[which(lassosol$cvm>cvtop)]
  #lambda_opt <- lambda_opt[1]
  #lambda_opt <- lambda_opt[1]
  #sol <- glmnet(W,-w,lambda=lambda_opt)
  sol <- glmnet(W,-w,lambda=0)
  sol <- sol$beta
  
  solution <- as.matrix(sol)
  if(index[iteration]==1){
    v <- c( 1 ,solution[index[iteration]:length(solution)])
  } else{
    v <- c(solution[1:(index[iteration]-1)], 1 ,solution[index[iteration]:length(solution)])
  }
  #v <- rep(0,m)
  #v[index[i]] <- 1
  #v[-index[i]] <- sol
  
  # 3rd step: find the other communities using deflation
  A_def <- A - t(t(v)) %*% t(v)
  
  return(list(A=A_def,vector=v))
}
