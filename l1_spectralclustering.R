# code l1 spectral clustering
library(glmnet)

# rendre propre le code (Aperturbed, ind)
# regarder les sorties du lasso (trop restrictif? 10% du max-min)
# shrinker quand?
l1_spectralclustering <- function(A,k,index){
  
  # inputs : A, k, m, S1,...,Sk
  m <- ncol(A)
  
  # 1 com
  res <- find_com(A,k,m,index,1)
  com <- res[[2]]
  
  for (i in (2:k)){
    res <- find_com(res[[1]],k,m,index,i) 
    com <- cbind(com,res[[2]])
  }
}

find_com <- function(A,k,m,index,i){
  # 1st step: SVD
  svd <- eigen(A) 
  eigenvalues<- sort(svd$values,index.return=TRUE)
  eigenvectors<- svd$vectors[,eigenvalues$ix]
  U <- t(eigenvectors[,1:(n-k+i-1)])
  
  # 2nd step: solve the lasso problem to find the first community
  w <- U[,index[i]]
  W <- U[,-index[i]]
  
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
  
  solution=as.matrix(sol)
  if(i==1){v=c( 1 ,solution[index[i]:length(solution)])}
  else{
    v=c(solution[1:(index[i]-1)], 1 ,solution[index[i]:length(solution)])
  }
  #v <- rep(0,m)
  #v[index[i]] <- 1
  #v[-index[i]] <- sol
  
  # 3rd step: find the other communities using deflation
  A_def <- A - t(t(v)) %*% t(v)
  
  return(list(A=A_def,v=v))
}
