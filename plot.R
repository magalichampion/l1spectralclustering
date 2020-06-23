# plot some graphics (l1-spectral clustering)
path <- "/Users/mchampion/Desktop/l1_spectralclustering/l1spectralclustering/"
load(paste0(path,"data/EasySimu.Rdata"))
load(paste0(path,"results/Results_easySimu_TI.Rdata"))
resultsTI <- results # avec les bons indics
load(paste0(path,"results/Results_easySimu.Rdata"))
resultsNTI <- results

results <- resultsNTI

# plot the graph (from data)
lay <- layout.fruchterman.reingold(results$Structure$graph)
plot(results$Structure$graph,layout=lay,vertex.color=rep("lightcyan",ncol(data$A)))

# plot the comm
col <- rainbow(ncol(results$comm),s=0.4)
vertex_color <- rep(0,nrow(results$comm))
for (i in (1:ncol(results$comm))){
  I <- which(results$comm[,i]!=0)
  vertex_color[I] <- col[i]
  vertex_color[which(results$comm==1,arr.ind = TRUE)[i,1]] <- rainbow(ncol(results$comm))[i]
}

List <- list()
for (i in (1:ncol(results$comm))){
  CommunMember <- which(results$comm[,i]!=0)
  List <- c(List,list(CommunMember))
}
plot(graph,vertex.color=vertex_color,layout=lay,mark.groups = List,mark.border=NA)

CreateAdj <- function(comm){
  # comm : indicateur de comm
  A <- matrix(0,ncol=nrow(comm),nrow=nrow(comm))
  for (i in (1:ncol(comm))){
    I <- which(comm[,i]!=0)
    A[I,I]<- 1
  }
  A <- A-diag(diag(A))
}
