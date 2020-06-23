n_small<- c(10,20,50)
n_large <- c(100,200,500)
k_small <- c(2,3,4,5)
#k_large <- c(5,10,20)
p_int <- c(0.01,0.2,0.25,0.50)
p_ext <- p_int
#=permn(p_int)




for(n in 1:length(n_small) ){
  for(k in 1:length(k_small)){
    for(i in 1:length(p_int)){
      for (j in 1:length(p_ext)){
        p_in <- p_int[i]
        p_ext <- p[[j]]
        graphs <- list()
        data <- CreateDataSet(k=k_small[k],n=n_small[n],p=list(p_inside=p_int[i],p_outside=p_ext[j]),print.plot = TRUE)
        ClustersLength <- data$ClustersLength
        graphs <-  c(graphs,list(data))
        for(c in 2:100){
          data <- CreateDataSet(k=k_small[k],n=n_small[n],p=list(p_inside=p_int[i],p_outside=p_ext[j]),print.plot = FALSE,ClustersLength = ClustersLength)
          
          graphs <-  c(graphs,list(data))
        }
        names(graphs) <- paste0("Graph",1:100)
        save(graphs,file=paste0("/Users/camille/Desktop/Algorithmes/Spectral_Clustering/data/data_n=",n_small[n],"_k=",k_small[k],"_p_inside=",p_int[i],"_p_outside=",p_ext[j],".Rdata"))
      }
    }
  }
}


T1<-Sys.time()
results <- l1_spectralclustering(M=data$A_hat)
T2<-Sys.time()

Tdiff= difftime(time2, time1)

######## attention afficher en plus de A clusterslength 
# refaire tourner le code avec 

ideal_com <- function(ClustersLength=data$ClustersLength,A=data$A){
  mat=data$A+diag(dim(data$A)[2])
  perfect <- matrix(0,sum(data$ClustersLength),length(ClustersLength))
  for(i in 1:length(data$ClustersLength)){
    perfect[,i]<- mat[,cumsum(data$ClustersLength)[i]]
  }
  return(perfect)
}
######ROC
library(ROCR)
real <- ideal_com[data$ClustersLength,data$A]
pred <- prediction(results,real)
pred

#######score 
scoring_detection=function(predicted,cluster){
  lclust <-length(cluster) 
  z <- length(which(predicted==0)) 
  if(z!=0){ 
    core <- which(predicted==0)
    cluster <- cluster[-core]
    predicted <- predicted[-core]
  }
  S <- matrix(0,ncol=length(unique(predicted)),nrow=1)
  for (i in 1:length(unique(predicted))){
    L <- which(predicted==i) 
    nbclust <- unique(cluster[L]) 
    if (length(nbclust)==1){
      S[i]=length(L)
    }
    else{
      l=matrix(0,ncol=length(nbclust),nrow=1)
      for (j in 1:length(nbclust)){
        l[j]=length(which(cluster[L]==nbclust[j])) 
      }
      max <- max(l) 
      S[i]=max
    }
  }
  Scoring <- (1/lclust)*sum(S)
  return(Scoring)
}

####temps de calcul
ptm=proc.time()
##algo a enregistrer
proc.time()-ptm

######spectral reality
spec_result<- scoring_detection(spec,s1)


######## core reality
s1 <- rep(1,83)
for (i in 2:100){
  s <- rep(i,83)
  s_tot <- c(s1,s)
  s1=s_tot
}
scoring_detection(core[,2],s1)
