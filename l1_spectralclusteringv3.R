library(glmnet)
library(igraph)
library(NMI)
library(purrr)
library(pracma)

# code l1 spectral clustering - v3
data <- CreateDataSet(k=3,n=10,p=list(p_inside=0.1,p_outside=0.1),print.plot = TRUE)
results <- l1_spectralclustering(M=data$A_hat, pen="lasso")

##### Create data set #####
CreateDataSet <- function(k,n,p,print.plot=TRUE,ClustersLength = NULL){
  # k: number of clusters
  # n: number of nodes
  # perturbation: vector of perturbations, p1 (inside) and p2 (outside)
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

##### l1-spectral clustering #####
l1_spectralclustering <- function(M, k = NULL, index = NULL, pen){
  # M: the matrix we are working on (adjacency matrix)
  # k: true number of clusters (not necessary needed)
  # index: nodes belonging to the communities (not necessary needed)
  
  # 1st step: find the components
  Structure <- FindStructure(M)
  
  # 2nd step: find the optimal number of clusters
  clusters <- FindNbrClusters(M,Structure  = Structure,k = k)
  
  # 3rd step: find the indices of the community
  Elements <- FindElement(M = M, Structure = Structure, k = clusters,index = index)
    
  # 4th step: the l1-spectral clustering algorithm (each component are treated independtly)
  comm <- matrix(0,nrow=ncol(M),ncol=clusters$nbr_clusters_total)
  S <- cumsum(unlist(clusters$nbr_clusters))
  for (i in (1:length(Structure$groups))){
    Mtmp <- M[Structure$groups[[i]],Structure$groups[[i]]]
    clusters_tmp <- clusters$nbr_clusters[[i]]
    indices_tmp <- Elements$indices[which(Elements$indices%in%Structure$groups[[i]])]
    indices_tmp <- match(indices_tmp,Structure$groups[[i]])
    score_tmp <- Elements$score[[i]]
    names(score_tmp) <- paste0("Node",match(as.numeric(substring(names(Elements$score[[i]]),5)),Structure$groups[[i]]))
    Elements_tmp <- list(score = score_tmp,indices = indices_tmp)
    
    results <- l1spectral(M = Mtmp, k = clusters_tmp, elements = Elements_tmp,pen=pen)
    
    if (!is.null(ncol(results))){
      if (ncol(results)!=clusters$nbr_clusters[[i]]){
        cluster <- clusters$nbr_clusters
        cluster[[i]] <- ncol(results)
        S <- cumsum(unlist(cluster))
      }
    }
    if (i==1){
      comm[Structure$groups[[i]],1:S[1]] <- results
    } else {
      comm[Structure$groups[[i]],(S[i-1]+1):S[i]] <- results
    }
  }
  if (length(which(colSums(comm)==0))>0){
    comm <- comm[,-which(colSums(comm)==0)]
  }
  if (!is.null(ncol(comm))){
    if (length(which(rowSums(comm>0)>1))){
      I <- which(rowSums(comm>0)>1)
      for (i in (1:length(which(rowSums(comm>0)>1)))){
        comm[I[i],-which.max(comm[I[i],])] <- 0
      }
    }
  }
  comm[comm>0] <- 1
  return(list(comm=comm,Structure=Structure,clusters=clusters,Elements =Elements))
}

l1spectral <- function(M, k, elements,pen){
  if (length(M)==1){
    # only one node in the community
    comm <- 1
  } else {
    # code for running the l1-spectral algorithm for one component
    indices <- elements$indices
  
    # 1st step: svd on M
    n <- ncol(M)
    svd <- eigen(M) 
    eigenvalues <- sort(svd$values,index.return=TRUE)
    eigenvectors <- svd$vectors[,eigenvalues$ix]
  
    # 2nd step: loop on the number of clusters  
    algo <- "stop"
    comm <- c()
    DoubleNodes <- c()
  
    while (algo == "stop"){
      if (length(DoubleNodes)>0){
        # find other indices
        I <- elements$score[-which(names(elements$score)%in%DoubleNodes)]
        if (length(I)==0){
          print("One cluster disappears.")
          algo <- "continue"
          break
        } else if (length(I)<k){
          k <- length(I)
        } else {
          I <- names(sort(I[1:k]))
          indices <- as.numeric(substring(I, 5))
        }
      }
    
      if (k>1){
        eigenvectors_tmp <- eigenvectors
        comm <- c()
        for (i in (1:k)){
          # 3rd step: check the indices (only if i>1)
          if (i>1){
            if (length(which(v[indices[-(i-1)]]>0))>0){
              print("Find other community indices.")
              doubleNodes <- paste0("Node",indices[-(i-1)][which(v[indices[-(i-1)]]>0)])
              DoubleNodes <- c(DoubleNodes,doubleNodes)
              algo <- "stop"
              break
            } else {
              algo <- "continue"
            }
          } 
    
          # 4th step: Gram-Schmidt (only if i>1)
          if (i>1){
            eigenvectors_tmp <- eigenvectors_tmp[,-(n-k+i-1)]
            eigenvectors_tmp <- cbind(v,eigenvectors_tmp)
      
            eigenvectors_tmp <- GramSchmidt(eigenvectors_tmp)
            eigenvectors_tmp <- cbind(eigenvectors_tmp[,2:(n-k+i-1)],eigenvectors_tmp[,1],eigenvectors_tmp[,(n-k+i):n])
          }
      
          # 5th step: solve the lasso
          U <- t(eigenvectors_tmp[,1:(n-k+i-1)])
          v <- Lasso(U, n, indices, iteration = i,pen=pen,k)
          print(paste0("Iteration ",i," done."))
      
          # 6th step: save the community index
          comm <- cbind(comm,v)
          
          if (i==k){
            # check the indices for the last time
            if (length(which(v[indices[-i]]>0))>0){
              print("Find other community indices.")
              doubleNodes <- paste0("Node",indices[-i][which(v[indices[-i]]>0)])
              DoubleNodes <- c(DoubleNodes,doubleNodes)
              algo <- "stop"
              break
            } else {
              algo <- "continue"
            }
          } 
        }
      } else {
        comm <- rep(1,n)
        algo <- "continue"
      }
    }
  }
  return(comm)
}

Lasso <- function(U, n, indices, iteration,pen,k){
  w <- U[,indices[iteration]]
  W <- matrix(U[,-indices[iteration]],ncol=(ncol(U)-1),nrow=nrow(U))
  
  if (sum(w)==0 || length(w)==1){
    print("There is only one node in this cluster.")
    # w is constant - no more nodes in the cluster
    v <- rep(0,n)
    v[indices] <- 1
  } else if (length(which(w!=0))==1){
    sol <- glmnet(W,-w,lambda=0,lower.limits=0)
    sol <- sol$beta
    solution <- as.matrix(sol)
    
    solution[solution<0.5] <- 0
    
    if(indices[iteration]==1){
      v <- c( 1 ,solution[indices[iteration]:length(solution)])
    } else{
      v <- c(solution[1:(indices[iteration]-1)], 1 ,solution[indices[iteration]:length(solution)])
    }
    
  } else {
    if (pen == "lasso"){
      lassosol <- cv.glmnet(W,-w)
    
      cvtop<- min(lassosol$cvm)+5*(max(lassosol$cvm)-min(lassosol$cvm))/100
      plot(lassosol)
      abline(h=cvtop)
    
      error <- min(lassosol$cvm[lassosol$cvm>cvtop])
      lambda_opt <- lassosol$lambda[which(lassosol$cvm==error)]
    
      sol <- glmnet(W,-w,lambda=lambda_opt)
      sol <- sol$beta
      solution <- as.matrix(sol)
    
      I <- which(solution!=0)
      sol2 <- glmnet(W,-w,lambda=0)
      sol2 <- sol2$beta
      sol2 <- as.matrix(sol2)
      solution[I] <- sol2[I]
    } else if (pen=="elastic"){
      lassosol <- cv.glmnet(W,-w,alpha=0.5)
      
      cvtop<- min(lassosol$cvm)+5*(max(lassosol$cvm)-min(lassosol$cvm))/100
      plot(lassosol)
      abline(h=cvtop)
      
      error <- min(lassosol$cvm[lassosol$cvm>cvtop])
      lambda_opt <- lassosol$lambda[which(lassosol$cvm==error)]
      
      sol <- glmnet(W,-w,lambda=lambda_opt,alpha=0.5)
      sol <- sol$beta
      solution <- as.matrix(sol)
      
      I <- which(solution!=0)
      sol2 <- glmnet(W,-w,lambda=0)
      sol2 <- sol2$beta
      sol2 <- as.matrix(sol2)
      solution[I] <- sol2[I]
    } else {
      #lassosol <- cv.glmnet(W,-w)
      
   #   lambda_opt <- lassosol$lambda[which(lassosol$nzero<=(n-k))[1]]
  #    cvtop <- lassosol$cvm[which(lassosol$nzero==floor((n-1)/k))[1]]
  #    plot(lassosol)
  #    abline(h=cvtop)
      
   #   sol <- glmnet(W,-w,lambda=lambda_opt)
   #   sol <- glmnet(W,-w,lambda=0)
  #    sol <- sol$beta
  #    solution <- as.matrix(sol)
      
   #   I <- which(solution!=0)
      sol2 <- glmnet(W,-w,lambda=0)
      sol2 <- sol2$beta
      sol2 <- as.matrix(sol2)
      sol2[which(sol2<0.5)] <- 0
      solution <- sol2
     # solution[I] <- sol2[I]
    }
    if(indices[iteration]==1){
      v <- c( 1 ,solution[indices[iteration]:length(solution)])
    } else if (indices[iteration]==n){
      v <- c(solution[1:(indices[iteration]-1)], 1)
    } else {
      v <- c(solution[1:(indices[iteration]-1)], 1 ,solution[indices[iteration]:length(solution)])
    }
  }
  v[v<0] <- 0
  return(v)
}

FindElement <- function(M, Structure, k, index = NULL){
  n <- ncol(M)
  
  between <- function(comm){
    graph_small <- graph_from_adjacency_matrix(M[comm,comm],mode="undirected")
    b <- betweenness(graph_small,directed = FALSE,normalized=TRUE)
    names(b) <- paste0("Node",comm)
    I <- order(b,decreasing = TRUE)
    b <- b[I]
    return(b)
  }
  betweenness <- lapply(Structure$groups,between)
  
  if (!is.null(index)){
    if (length(index)!=k$nbr_clusters_total){
      print(paste0("The number of communities indices do not coincide with the number of communities. The algorithm will compute new indices for the ",clusters$nbr_clusters_total," communities."))
      continue <- TRUE
    } else {
      continue <- FALSE
      indices <- index
    }
  } else {
    continue <- TRUE
  }
  if (continue == TRUE){
    if (length(betweenness)!=1 && length(betweenness) != (length(k$nbr_clusters)-1)){
      print(paste0("Please, be careful. There is a non-negative number of connected components in the graph and you set the number of clusters yourself. The choice of indices may be wrong as the algorithm will randomly choose the same number of indices for each community."))
      
      d <- k$nbr_clusters_total%/%length(Structure$groups)
      r <- k$nbr_clusters_total%%length(Structure$groups)
      K <- rep(d,length(Structure$groups))
      K[(length(K)-r+1):length(K)] <- K[(length(K)-r+1):length(K)]+1
      
      Nodes <- c()
      for (i in (length(betweenness):1)){ 
        b <- betweenness[[i]][1:K[i]]
        Nodes <- c(names(rev(b)),Nodes)
      }
      Nodes <- as.numeric(substring(Nodes, 5))
    } else {
      Nodes <- c()
      for (i in (1:length(betweenness))){ 
        b <- betweenness[[i]][1:unlist(k$nbr_clusters)[i]]
        Nodes <- c(Nodes,names(rev(b)))
      }
      Nodes <- as.numeric(substring(Nodes, 5))
    }
    
    print(paste0(c("The",k$nbr_clusters_total,"indices of the communities are",Nodes,"."), collapse=" "))
  } else {
    Nodes <- index
    print(paste0(c("The (provided)",k$nbr_clusters_total,"indices of the communities are",Nodes,"."), collapse=" "))
  }
  return(list(score=betweenness,indices=Nodes))
}

FindNbrClusters <- function(M,Structure, k = NULL){
  nbr_group <- length(Structure$groups)
 
  Eigen_list <- function(group){
    if (length(group)>1){
      D <- diag(rowSums(M[group,group])) # degree matrix
      L <- D-M[group,group]
      svd <- eigen(L) 
      eigenvalues <- sort(svd$values)
    } else {
      eigenvalues <- NA
    }
    return(eigenvalues)
  }
  
  Gap <- function(eigenvalue){
    if (length(eigenvalue)==1){
      nbr_cluster <- 1
    } else {
      gap <- c()
      for (i in (2:length(eigenvalue))){
        gap <- c(gap,eigenvalue[i]-eigenvalue[i-1])
      }
      gap_ecart <- c(gap,0)-c(0,gap)
      nbr_cluster <- which(gap_ecart[2:length(gap_ecart)]>0.20)[1]+1
      if (is.na(nbr_cluster)){
        nbr_cluster=1
      }
    }
    return(nbr_cluster)
  }
  
  if (is.null(k)){
    if (nbr_group>1){
      eigenvalues <- lapply(X = c(Structure$groups,list(all=c(1:ncol(M)))),FUN = Eigen_list)
    } else {
      eigenvalues <- list(all=Eigen_list(Structure$groups[[1]]))
    } 
    par(mfrow=c(1,1))
    plot(eigenvalues$all,main="Eigenvalues of the Laplacian matrix",ylab="Eigenvalues",xlab="",type="b")
    if (nbr_group>1){
      for (i in (1:nbr_group)){
        points(eigenvalues[[i]],col=rainbow(nbr_group)[i],type="b")
      }
      legend("bottomright",legend = c("All nodes",paste0("Connected component ",c(1:nbr_group))),col=c("black",rainbow(nbr_group)[1:nbr_group]),lty=1)
    }

    gaps <- lapply(X = eigenvalues,FUN = Gap)
    par(xpd=FALSE) 
    for (i in (1:length(gaps))){
      color <- c(rainbow(nbr_group),"black")
      abline(v=gaps[[i]],col=color[i],ylim=c(0,10),lty=3)
    }
    
    if (length(gaps)==1){
      nbr_clusters_total <- gaps$all
      print(paste0("The optimal number of clusters is ",gaps$all,"."))
    } else {
      nbr_clusters_total <- sum(unlist(gaps)[1:(length(gaps)-1)])
      print(paste0("The optimal number of clusters is ",nbr_clusters_total,"."))
      print(paste0(c("Here,",nbr_group,"connected components were detected. Each of them should be clustered into",unlist(gaps[1:nbr_group]),"clusters."),collapse=" "))
    }
    
  } else if (nbr_group==1 && !is.null(k)){
    nbr_clusters_total <- k
    gaps <- list(all=k)
    print(paste0("The provided number of clusters is ",nbr_clusters_total,"."))
    
  } else {
    # k is provided but there are more than one component
    eigenvalues <- lapply(X = c(Structure$groups,list(all=c(1:ncol(M)))),FUN = Eigen_list)
    par(mfrow=c(1,1))
    plot(eigenvalues$all,main="Eigenvalues of the Laplacian matrix",ylab="Eigenvalues",xlab="",type="b")
    for (i in (1:nbr_group)){
      points(eigenvalues[[i]],col=rainbow(nbr_group)[i],type="b")
    }
    legend("bottomright",legend = c("All nodes",paste0("Connected component ",c(1:nbr_group))),col=c("black",rainbow(nbr_group)[1:nbr_group]),lty=1)

    gaps <- lapply(X = eigenvalues,FUN = Gap)
    par(xpd=FALSE) 
    for (i in (1:length(gaps))){
      color <- c(rainbow(nbr_group),"black")
      abline(v=gaps[[i]],col=color[i],ylim=c(0,10),lty=3)
    }
    
    if (sum(unlist(gaps)[1:(length(gaps)-1)]) != k){
      # there is a problem
      gaps[[(length(gaps)-1)]] <- k-sum(unlist(gaps)[1:(length(gaps)-2)])
    }
    nbr_clusters_total <- sum(unlist(gaps)[1:(length(gaps)-1)])
    print(paste0("The provided number of clusters is ",k,"."))
    print(paste0(c("Here,",nbr_group,"connected components were detected. Each of them should be clustered into",unlist(gaps[1:nbr_group]),"clusters."),collapse=" "))
  }
  return(list(nbr_clusters=gaps,nbr_clusters_total=nbr_clusters_total))
}

FindStructure <- function(M){
  # function to explore the graph structure
  
  # define the graph
  graph <- graph_from_adjacency_matrix(M,mode="undirected")
  
  # find the connected components
  clu <- components(graph)
  groups <- groups(clu)
  if (length(groups)>1){
    groups <- groups[order(sapply(groups,length))]
  } else {
    groups <- groups
  }
  nbr_comp <- length(groups)
  names(groups) <- paste0("Component",c(1:nbr_comp))
  
  if (nbr_comp==1){
    print("The graph has only one connected component.")
  } else {
    print(paste0("The graph have ",nbr_comp," connected components. Each of them will be clustered using the l1-spectral clustering."))
  }
  return(list(graph=graph,groups=groups))
}

GramSchmidt <- function(U){
  #v <- U
  #B <- matrix(0,dim(U)[1],dim(U)[2])
  #B[,1]=v[,1]/sqrt(sum(v[,1]^2))
  #for(i in 2:(dim(U)[2])){
  #  res=v[,i]
    
  #  for(j in 1:(i-1)){
  #    res=res-(crossprod(B[,j],v[,i])*B[,j]/sqrt(sum(B[,j]^2)))
  #  }
    
  #  B[,i]=res
  #  B[,i]=B[,i]/sqrt(sum(B[,i]^2))
  #}
  #B <- gramSchmidt(U)
  B <- grahm_schimdtCpp(A)
  B <- B$Q
  return(B)
}
