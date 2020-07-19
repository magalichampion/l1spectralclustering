n_small<- c(10,20,50,100)
k_small <- c(2,3,4,5,10)
p_int <- c(0.01,0.1,0.2,0.25,0.5)
p_ext <- p_int

path <- "/Users/mchampion/Desktop/l1_spectralclustering/l1spectralclustering/"
# path <- "/Users/camille/Desktop/Algorithmes/Spectral_Clustering/"

# creation des graphes
for(nit in 1:length(n_small)){
  for(kit in 1:length(k_small)){
    for(iit in 1:length(p_int)){
      for (jit in 1:length(p_ext)){
        graphs <- list()
        data <- CreateDataSet(k=k_small[kit],n=n_small[nit],p=list(p_inside=p_int[iit],p_outside=p_ext[jit]),print.plot = TRUE)
        ClustersLength <- data$ClustersLength
        graphs <-  c(graphs,list(data))
        for(c in 2:100){
          data <- CreateDataSet(k=k_small[kit],n=n_small[nit],p=list(p_inside=p_int[iit],p_outside=p_ext[jit]),print.plot = FALSE,ClustersLength = ClustersLength)
          
          graphs <-  c(graphs,list(data))
        }
        names(graphs) <- paste0("Graph",1:100)
        save(graphs,file=paste0(path,"data/data_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
      }
    }
  }
}

# run the code
for(nit in 1:length(n_small) ){
  for(kit in 1:length(k_small)){
    for(iit in 1:length(p_int)){
      for (jit in 1:length(p_ext)){
        load(paste0(path,"data/data_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
        Results <- list()
        Tdiff <- c()
        NMI <- c()
        for (l in 1:100){
          graph <- graphs[[l]]
          T1 <- Sys.time()
          results <- l1_spectralclustering(M=graph$A_hat,pen="autre")
          T2 <- Sys.time()
          tdiff= difftime(T2, T1)
          Tdiff <- c(Tdiff,tdiff)
          nmi <- Compute_NMIscore(comm=results$comm,A=graph$A,method="l1")
          NMI <- c(NMI,nmi)
          Results <- c(Results,list(results))
        }
        Results <- c(Results,list(Tdiff),list(NMI))
        names(Results) <- c(paste0("Graphs",1:100),"ProcTime","NMI")
        save(Results,file=paste0(path,"results/Results_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
      }
    }
  }
}

# l1 - with predefined indices
for(nit in 1:length(n_small) ){
  for(kit in 1:length(k_small)){
    for(iit in 1:length(p_int)){
      for (jit in 1:length(p_ext)){
        load(paste0(path,"data/data_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
        Results <- list()
        Tdiff <- c()
        NMI <- c()
        for (l in 1:100){
          graph <- graphs[[l]]
          UneBoucle <- function(graph,K=NULL){
            ScoreNMI <- c()
            Tps_total <- c()
            Community_total <- list()
          
            StructureReal <- FindStructure(graph$A)
            Structure <- FindStructure(graph$A_hat)
          
            # if more than one component, define a number of clusters by component
            Comm <- list()
            Tps <- list()
            if (is.null(K)){
              K <- c()
            }
            for (i in (1:length(Structure$groups))){
              indices <- Structure$groups[[i]]
              
              if (length(indices)>1){
                loc <- c()
                for (j in (1:length(indices))){
                  I <- sapply(StructureReal$groups, FUN = function(x){
                    I <- length(which(x==indices[[j]]))
                  })
                  loc <- c(loc,which(I>0))
                }
                if (is.null(K)){
                  k <- length(unique(loc))
                  K <- c(K,k)
                } else {
                  k <- K[i]
                }
            
                Mtmp <- graph$A_hat[Structure$groups[[i]],Structure$groups[[i]]]
                Mtmp_real <- graph$A[Structure$groups[[i]],Structure$groups[[i]]]
                StructureReal_tmp <- FindStructure(Mtmp_real)
            
                comb <- StructureReal_tmp$groups %>%
                    cross()
                indicesBis <- list()
                Rand <- sample(1:length(comb),100)
                #for (t in (1:length(comb))){
                for (u in (1:length(Rand))){
                  t <- Rand[u]
                  indicesBis <- c(indicesBis,list(unlist(comb[[t]])))
                }
            
                results_tmp <- c()
                tdiff_tmp <- c()
                nmi_tmp <- c()
                
                for (t in (1:length(indicesBis))){
                  T1 <- Sys.time()
                  results <- l1_spectralclustering(M=Mtmp,pen="autre",k = k,index = indicesBis[[t]])
                  if (ncol(results$comm)!=k){
                    comm_tmp <- matrix(0,nrow=nrow(results$comm),ncol=k)
                    comm_tmp[,1:ncol(results$comm)] <- results$comm
                    results$comm <- comm_tmp
                  }
                  T2 <- Sys.time()
                  tdiff = difftime(T2, T1)
                  tdiff_tmp <- c(tdiff_tmp,tdiff)
                  results_tmp <- c(results_tmp,list(results$comm))
                } 
                Tps <- c(Tps, list(tdiff_tmp))
                Comm <- c(Comm,list(results_tmp))
              } else {
                K <- c(K,1)
                indicesBis <- list(1)
                Comm <- c(Comm,list(rep(1,1)))
                Tps <- c(Tps,list(0))
              }
            }
            names(Tps) = names(Comm) <- names(Structure$groups)
            
            if (length(Comm)>1){
              comb <- Structure$groups %>%
                cross()
              indicesBis <- list()
              for (t in (1:length(comb))){
                indicesBis <- c(indicesBis,list(unlist(comb[[t]])))
              }
          
              for (t in (1:length(indicesBis))){
                Community <- matrix(0,nrow=ncol(graph$A),ncol=sum(K))
                S <- indicesBis[[t]]
                tps_tmp <- 0
                for (s in (1:length(S))){
                  I <- match(S[s],Structure$groups[[s]])
                  res <- Comm[[s]][[I]]
                  if (!is.null(ncol(res))){
                    if (ncol(res)!=K[s]){
                      res <- cbind(res,rep(0,nrow(res)))
                    }
                    if (s==1){
                      Community[Structure$groups[[s]],1:ncol(res)] <- res
                    } else {
                      Community[Structure$groups[[s]],(cumsum(K)[s-1]+1):cumsum(K)[s]] <- res
                    }
                  } else {
                    if (s==1){
                      Community[Structure$groups[[s]],1] <- res
                    } else {
                      Community[Structure$groups[[s]],(cumsum(K)[s-1]+1):cumsum(K)[s]] <- res
                    }
                  }
                  tps_tmp <- tps_tmp + Tps[[s]][[I]]
                }
                score <- Compute_NMIscore(comm=Community,A=graph$A,method="l1")
                ScoreNMI <- c(ScoreNMI,score)
                Tps_total <- c(Tps_total,tps_tmp)
                Community_total <- c(Community_total,list(Community))
              }
            } else {
              for (s in (1:length(Comm[[1]]))){
                Community <- Comm[[1]][[s]]
                score <- Compute_NMIscore(comm=Community,A=graph$A,method="l1")
                ScoreNMI <- c(ScoreNMI,score)
                Tps_total <- c(Tps_total, Tps[[1]][[s]])
                Community_total <- c(Community_total,list(Community))
              }
            }
            I <- which.max(ScoreNMI)
            return(list(Results=Community_total[[I]],Tdiff=Tps_total[I],NMI=ScoreNMI[I]))
          }
          UneBoucle2 <- possibly(UneBoucle, otherwise = list(Results=NA,Tdiff=NA,NMI=NA))
          
          results <- UneBoucle2(graph)
          Results <- c(Results,list(results$Results)) 
          Tdiff <- c(Tdiff,results$Tdiff)
          NMI <- c(NMI,results$NMI)
        }
        Results <- c(Results,list(Tdiff),list(NMI))
        names(Results) <- c(paste0("Graphs",1:100),"ProcTime","NMI")
        save(Results,file=paste0(path,"results/IndicResults_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
      }
    }
  }
}
        
# spectral clustering
library(anocva)
for(nit in 1:length(n_small) ){
  for(kit in 1:length(k_small)){
    for(iit in 1:length(p_int)){
      for (jit in 1:length(p_ext)){
        load(paste0(path,"data/data_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
        Results <- list()
        Tdiff <- c()
        NMI <- c()
        for (l in (1:100)){
          graph <- graphs[[l]]
          T1 <- Sys.time()
          spectral <- spectralClustering(graph$A_hat,k=length(graph$ClustersLength))
          T2 <- Sys.time()
          tdiff = difftime(T2, T1)
          Tdiff <- c(Tdiff,tdiff)
          score <- Compute_NMIscore(comm=spectral,A=graph$A,method="spectral")
          NMI <- c(NMI,score)
          Results <- c(Results,list(spectral))
        }
        Results <- c(Results,list(Tdiff),list(NMI))
        names(Results) <- c(paste0("Graphs",1:100),"ProcTime","NMI")
        save(Results,file=paste0(path,"results/SpectralClustering_Results_n=",n_small[nit],"_k=",k_small[kit],"_p_inside=",p_int[iit],"_p_outside=",p_ext[jit],".Rdata"))
      }
    }
  }
}

# recuperer la taille des clusters
load(paste0(path,"results/Results_n=",10,"_k=",2,"_p_inside=",0.1,"_p_outside=",0.1,".Rdata"))
k <- c()
for (l in (1:100)){
  k <- c(k,ncol(Results[[l]]$comm))
}

load("~/Desktop/l1_spectralclustering/l1spectralclustering/results/Scores.rdata")
scores <- c(scores[,1],scores[,2],score2)
scores <- data.frame(cbind(scores=scores,method=c(rep("l1",100),rep("spectral",100),rep("indic",100)),clusters=c(k,rep(2,200))))
scores$scores <- as.numeric(as.character(scores$scores))
scores$clusters <- as.numeric(as.character(scores$clusters))
boxplot(scores~method,data=scores)
I <- which(scores$clusters==2)      
boxplot(scores~method,data=scores[I,])

load("~/Dropbox/L1Spectral/results/Results_n=50_k=5_p_inside=0.25_p_outside=0.25.Rdata")
ResultsL1 <- Results$NMI
k <- c()
for (l in (1:100)){
  k <- c(k,ncol(Results[[l]]$comm))
}
load("~/Dropbox/L1Spectral/results/SpectralClustering_Results_n=50_k=5_p_inside=0.25_p_outside=0.25.Rdata")
ResultsSpec <- Results$NMI

scores <- data.frame(cbind(scores=c(ResultsL1,ResultsSpec),method=c(rep("l1",100),rep("spectral",100)),clusters=c(k,rep(5,100))))
scores$scores <- as.numeric(as.character(scores$scores))
scores$clusters <- as.numeric(as.character(scores$clusters))
boxplot(scores~method,data=scores)

I <- which(scores$clusters==5)      
boxplot(scores~method,data=scores[I,])
