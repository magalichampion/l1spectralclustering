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

recovery_adjacency <- function(A,com){
  
  li1<- c()
  m=dim(com)[1]
  for (i in 1:k){
    ind1 <- which(com[,i]==1)
    li1[i] <- length(ind1)
  }
  adja<- matrix(0,m,m)

  mat0 <- matrix(c(rep(com[,1],li1[1])),li1[1],dim(com)[1],byrow = TRUE)
  for(j in 2:(dim(com)[2])){
    mat <- matrix(c(rep(com[,j],li1[j])),li1[j],dim(com)[1],byrow = TRUE)
    mat2 <- rbind(mat0,mat)
    mat0=mat2
  }
  mat0[lower.tri(mat0)] <- 0
  mat0=mat0+t(mat0)
  mat0=mat0-diag(diag(mat0))
   # i1<- c()
   # m=dim(com)[1]
   # ind0 <- which(com[,1]==1)
   # for (i in 2:k){
   #   ind <- which(com[,i]==1)
   #   ind0 <- c(ind0,ind)
   # }
   # uni<- unique(ind0)
   # 
   # #indices0 <- which(ind0==unique(ind0)[1])
   # for (j in 1:length(uni)){
   #   #indices <- which(ind0==unique(ind0)[j]) 
   #   #indices1 <- c(indices0,indices)
   #   #indices0 <-indices1 
   #   indices0[j] <- which(ind0==unique(ind0)[j])[1]
   # }
   #   inl <- sort(indices0,index.return=TRUE)
   #   uni <- uni[inl$ix]
   #   A_recovery  <- A[uni,uni]
  return(A_recovery)
}