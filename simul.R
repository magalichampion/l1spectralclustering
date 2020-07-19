
files <- list.files(here("results"))
nb_files <- length(files_names)
data_names <- vector("list",length=nb_files)

files_names <- list.files(path=here("results"))
files_names
nb_files <- length(files_names)
#data_names <- vector("list",length=nb_files)

files_names2 <- list.files(path=here("data"))
files_names2
nb_files2 <- length(files_names2)

for (i in 1:nb_files){
  load(files_names[i])
  real<-ideal_com(ClustersLength=data$ClustersLength,A=data$A)
  pred <- prediction(results,real)
  clusters_pred <- clust_appt(ideal=results)
  clusters_real <- clust_appt(ideal=real)
  score <- NMI(cbind(c(1:length(clusters_pred)),clusters_pred),cbind(c(1:length(clusters_real)),clusters_real))
}

load("Results_n=10_k=2_p_inside=0.1_p_outside=0.1.Rdata")
Results_l1 <- Results
load("data_n=10_k=2_p_inside=0.1_p_outside=0.1.Rdata")
load("SpectralClustering_Results_n=10_k=2_p_inside=0.1_p_outside=0.1.Rdata")
Results_spec <- Results

n_small<- c(10,20,50)
n_large <- c(100,200,500)
k_small <- c(2,3,4,5)
k_large <- c(5,10,20)
p_int <- c(0.01,0.1,0.2,0.25,0.5)
p_ext <- p_int

for(n in 1:length(n_small) ){
  for(k in 1:length(k_small)){
    for(i in 1:length(p_int)){
      for (j in 1:length(p_ext)){
        # data
        load(paste0("/Users/camille/Documents/articlespectral/l1spectralclustering/data/data_n=",n_small[n],"_k=",k_small[k],"_p_inside=",p_int[i],"_p_outside=",p_ext[j],".Rdata"))
        
        # spectral
        load(paste0("/Users/camille/Documents/articlespectral/l1spectralclustering/results/SpectralClustering_Results_n=",n_small[n],"_k=",k_small[k],"_p_inside=",p_int[i],"_p_outside=",p_ext[j],".Rdata"))
        Results_spec <- Results
        
        # l1
        load(paste0("/Users/camille//Documents/articlespectral/l1spectralclustering/results/Results_n=",n_small[n],"_k=",k_small[k],"_p_inside=",p_int[i],"_p_outside=",p_ext[j],".Rdata"))
        Results_l1 <- Results
        
        # real
        real<-ideal_com(ClustersLength=graphs[[1]]$ClustersLength,A=graphs[[1]]$A)
        clusters_real <- clust_appt(ideal=real)
        
        # score
        score<- c()
        score_spec <- c()
        for(l in 1:100){
          if(is.na(Results_l1[[l]])==TRUE){
            score<-c(score,NA)
          }else{
          clusters_pred <- clust_appt(ideal=Results_l1[[l]]$comm)
          score<- c(score,NMI(cbind(c(1:length(clusters_pred)),clusters_pred),cbind(c(1:length(clusters_real)),clusters_real))$value)
          score_spec <- c(score_spec,NMI(cbind(c(1:length(Results_spec[[l]])),Results_spec[[l]]),cbind(c(1:length(clusters_real)),clusters_real))$value)
          }
          }

        Results_l1<- c(Results_l1,list(NMI=score))
        Results_spec<- c(Results_spec,list(NMI=score_spec))

        Results <- Results_l1
        save(Results,file=paste0("/Users/camille//Documents/articlespectral/l1spectralclustering/results/Results_n=",n_small[n],"_k=",k_small[k],"_p_inside=",p_int[i],"_p_outside=",p_ext[j],".Rdata"))
        
        Results <- Results_spec
        save(Results,file=paste0("/Users/camille//Documents/articlespectral/l1spectralclustering/results/SpectralClustering_Results_n=",n_small[n],"_k=",k_small[k],"_p_inside=",p_int[i],"_p_outside=",p_ext[j],".Rdata"))
        
      }
    }
  }
}


real<-ideal_com(ClustersLength=graphs[[1]]$ClustersLength,A=graphs[[1]]$A)
clusters_real <- clust_appt(ideal=real)
score<- c()
score_spec <- c()
for(i in 1:100){
   clusters_pred <- clust_appt(ideal=Results_l1[[i]]$comm)
   score[i] <- NMI(cbind(c(1:length(clusters_pred)),clusters_pred),cbind(c(1:length(clusters_real)),clusters_real))
   score_spec[i] <- NMI(cbind(c(1:length(Results_spec[[i]])),Results_spec[[i]]),cbind(c(1:length(clusters_real)),clusters_real))
}
score_l1=unlist(score)
score_spectral <- unlist(score_spec)
scores<-cbind(score_l1,score_spectral)
save(scores,file="Scores.rdata")
