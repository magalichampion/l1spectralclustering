merge_simulation_results <- function(path, n = NULL, d = NULL, methods=NULL) {
  all_files <- list.files(paste0(path,"results"), full.names = FALSE)
  
  Data_Full <- data.frame()
  Data_Summary <- data.frame()

  for (p_inside in p_inside_range){
    for (p_outside in p_outside_range){
      data_pattern <- paste0(path,"data/data_n=", n,"_density=",sprintf("%.2f", d), "_p_inside=", sprintf("%.2f", p_inside), "_p_outside=", sprintf("%.2f", p_outside),".Rdata")
      load(data_pattern)
      
      for (current_method in methods){
        file_pattern <- paste0("results_", current_method, 
                               "_n=", n, 
                               "_density=", sprintf("%.2f", d), 
                               "_p_inside=", sprintf("%.2f", p_inside), 
                               "_p_outside=", sprintf("%.2f", p_outside))
        
        target_file <- all_files[str_detect(all_files, file_pattern)]
        
        if (length(target_file) > 0) {
          load(file.path(paste0(path, "results/",target_file[1])))
          
          AMI <- results$AMI
          
          miss_counts <- c()
          AMI_corrected <- c()
          coverage <- c()
          for (g in (1:(length(results)-4))) {
            res <- results[[g]]
            
            clus_est <- res$clusters
            
            if (!any(is.na(clus_est))){
              true_clus <- igraph::components(graph_from_adjacency_matrix(graphs[[g]]$A, mode="undirected"))$membership
              
              idx <- which(clus_est != 0) # Find nodes your method actually clustered
              if(length(idx) > 0) {
                AMI_corrected <- c(AMI_corrected,aricode::AMI(clus_est[idx], true_clus[idx]))
                coverage <- c(coverage,length(idx) / length(true_clus))
              } else {
                AMI_corrected <- c(AMI_corrected,0)
                coverage <- c(coverage,0)
              }
              
              t <- table(as.vector(clus_est), true_clus)
              if (rownames(t)[1]=="0"){
                t <- t[-1,]
              }
              # include the number of missclassified nodes here!!!
              errors <- 0
            
              if (nrow(t) > 0) {
                best_match <- apply(t, 1, which.max)
              
                correct_nodes <- sum(sapply(1:nrow(t), function(row_idx) {
                  t[row_idx, best_match[row_idx]]
                }))
              
                current_miss <- errors + (sum(t) - correct_nodes)
              }
            } else {
              current_miss <- NA
              AMI_corrected <- c(AMI_corrected,NA)
              coverage <- c(coverage,NA)
            }
            miss_counts <- c(miss_counts, current_miss)
          }
          
          Data_Full <- rbind(Data_Full, data.frame(
            AMI    = AMI,
            AMI_corrected = AMI_corrected,
            coverage = coverage,
            Miss = miss_counts,
            Miss_perc = miss_counts / n,
            Method = current_method,
            p_inside   = p_inside,
            p_outside  = p_outside,
            n      = n,
            density = d
          ))
          
          Data_Summary <- rbind(Data_Summary, data.frame(
            Method     = current_method,
            Median_AMI = median(AMI, na.rm = TRUE),
            Mean_AMI   = mean(AMI, na.rm = TRUE),
            Mean_Miss = mean(miss_counts, na.rm = TRUE),
            Mean_Miss_perc = mean(miss_counts, na.rm = TRUE) / n,
            p_inside   = p_inside,
            p_outside  = p_outside
          ))
        }
      }
    }
  }
  Data_Full <- Data_Full %>%
    mutate(Method = factor(Method, levels = c(
      "l1Spectral", "Spectral", "regSpectral", "robustSpectral", "GNN",
      "ST_l1Spectral", "ST_Spectral","Hybrid", "MCL"
    ))) %>%
    mutate(Category = case_when(
      Method %in% c("l1Spectral", "Spectral", "regSpectral", "robustSpectral","GNN") ~ "Non self-tuned",
      Method %in% c("ST_l1Spectral", "ST_Spectral","Hybrid", "MCL") ~ "Self-tuned",
    ))
  Data_Summary <- Data_Summary %>%
    mutate(Method = factor(Method, levels = c(
      "l1Spectral", "Spectral", "regSpectral", "robustSpectral", "GNN",
      "ST_l1Spectral", "ST_Spectral","Hybrid", "MCL"
    ))) %>%
    mutate(Category = case_when(
      Method %in% c("l1Spectral", "Spectral", "regSpectral", "robustSpectral","GNN") ~ "Non self-tuned",
      Method %in% c("ST_l1Spectral", "ST_Spectral","Hybrid", "MCL") ~ "Self-tuned",
    ))
  return(list(Full=Data_Full,Summary=Data_Summary))
}

plot_perf <- function(df, n, k) {
  # df is the output of merge_simulation_results function 
  p <- ggplot(df$Full,aes(x=Method, y=AMI,fill=Method)) +
    geom_boxplot(outlier.shape = NA)+scale_fill_manual(values=c("peachpuff","lightsalmon","lightsalmon3","tomato4","slategray2","steelblue1","steelblue3","steelblue4"))+
    facet_grid(p_outside~p_inside)+ scale_y_continuous(limits=c(0,1)) +
    theme_bw()+scale_x_discrete(guide = guide_axis(angle = 45))
  p <- p +geom_point(data = subset(df$Summary, Category == "Non self-tuned"),
                     mapping=aes(y=Mean_Miss_perc,x=Method),
                     stat="identity",
                     #position= position_dodge(width = .9),
                     pch=5,
                     size=1,color="red") + ggtitle(paste0("n=",n," and d=",d))
}

ggplot(subset(df$Full, Category == "Non self-tuned"), aes(x = coverage, y = AMI_corrected, color = Method)) +
  # Add a bit of jitter so overlapping points (like Spectral at 100%) are visible
  geom_jitter(alpha = 0.4, size = 1.5, width = 0.01, height = 0.01) +
  # Add a mean point for each method to highlight the center of the cluster
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, stroke = 1.5) +
  facet_grid(p_outside ~ p_inside, labeller = label_both) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
  theme_bw() +
  labs(
    title = "Precision-Recall Trade-off in Community Detection",
    subtitle = "Diamonds represent the centroid of each method's performance",
    x = "Coverage (% of Nodes Clustered)",
    y = "Precision (AMI on Clustered Nodes Only)"
  ) +
  scale_color_manual(values = c(
    "l1Spectral" = "red", "ST_l1Spectral" = "darkred",
    "Spectral" = "steelblue", "regSpectral" = "blue", 
    "GNN" = "darkgreen", "MCL" = "orange"
  ))