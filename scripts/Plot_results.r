#' Merge simulation results across parameter grids
#'
#' @param path Character. The base directory path containing "data/" and "results/" subdirectories.
#' @param n Integer. The number of nodes in the graph.
#' @param d Numeric. The density of the graph.
#' @param p_inside_range Numeric vector. Probabilities of internal edge perturbation.
#' @param p_outside_range Numeric vector. Probabilities of external edge perturbation.
#' @param method Character vector. Names of the clustering methods to aggregate.
#'
#' @return A list containing two data frames: 
#' \item{Full}{Every iteration result with metrics per graph.}
#' \item{Summary}{Mean/median metrics grouped by method and perturbation level.}
#' 
#' @export
#' @import dplyr
#' @import aricode
#' @import stringr
merge_simulation_results <- function(path, n = NULL, d = NULL, p_inside_range = c(0.01,0.1,0.25,0.5), p_outside_range = c(0.01,0.1,0.25,0.5), method = c("l1Spectral", "Spectral", "regSpectral","robustSpectral", "VGAE", "ST_l1Spectral", "Spectrum", "Hybrid", "MCL")) {
  all_files <- list.files(paste0(path,"/results"), full.names = FALSE)
  
  data_full_list <- list()

  for (p_inside in p_inside_range){
    for (p_outside in p_outside_range){
      data_pattern <- paste0(path,"/data/data_n=", n,"_density=",sprintf("%g", d), "_p_inside=", sprintf("%.2f", p_inside), "_p_outside=", sprintf("%.2f", p_outside),".Rdata")
      if (!file.exists(data_pattern)) next
      load(data_pattern)
      
      for (current_method in method){
        file_pattern <- paste0("results_", current_method, 
                               "_n=", n, 
                               "_density=", sprintf("%g", d), 
                               "_p_inside=", sprintf("%.2f", p_inside), 
                               "_p_outside=", sprintf("%.2f", p_outside))
        
        target_file <- all_files[str_detect(all_files, file_pattern)]
        
        if (length(target_file) > 0) {
          load(file.path(paste0(path, "/results/",target_file[1])))
          
          num_graphs <- length(results) - 4
          
          for (g in (1:num_graphs)) {
            res <- results[[g]]
            
            clus_est <- res$clusters
            true_clus <- graphs[[g]]$clusters
            
            if (any(is.na(clus_est))) {
              ami_corr <- NA
              cov_val <- NA
              current_miss <- NA
            } else {
              idx <- which(clus_est != 0) 
              if(length(idx) > 0) {
                ami_corr <- aricode::AMI(clus_est[idx], true_clus[idx])
                cov_val <- length(idx) / length(true_clus)
              } else {
                ami_corr <- 0
                cov_val <- 0
              }
              
              t <- table(as.vector(clus_est), true_clus)
              if (rownames(t)[1]=="0"){
                t <- t[-1,]
              }
              errors <- 0
            
              if (nrow(t) > 0) {
                best_match <- apply(t, 1, which.max)
              
                correct_nodes <- sum(sapply(1:nrow(t), function(row_idx) {
                  t[row_idx, best_match[row_idx]]
                }))
              
                current_miss <- errors + (sum(t) - correct_nodes)
              }
            }
          
            data_full_list[[length(data_full_list) + 1]] <- data.frame(
              AMI = results$AMI[g],
              AMI_corrected = ami_corr,
              coverage = cov_val,
              Miss = current_miss,
              Miss_perc = current_miss / n,
              Method = current_method,
              p_inside = p_inside,
              p_outside = p_outside,
              n = n,
              density = d
            )
          }
        }
      }
    }
  }
          
  Data_Full <- dplyr::bind_rows(data_full_list) %>%
      mutate(Method = factor(Method, levels = method)) %>%
      mutate(Category = case_when(
        Method %in% c("l1Spectral", "Spectral", "regSpectral", "robustSpectral", "VGAE") ~ "Non self-tuned",
          TRUE ~ "Self-tuned"
      ))
          
  Data_Summary <- Data_Full %>%
    group_by(Method, Category, p_inside, p_outside) %>%
      summarise(Median_AMI = median(AMI, na.rm = TRUE),
                Mean_AMI   = mean(AMI, na.rm = TRUE),
                Mean_Miss  = mean(Miss, na.rm = TRUE),
                Mean_Miss_perc = mean(Miss_perc, na.rm = TRUE),
                Mean_Coverage = mean(coverage, na.rm = TRUE),
                .groups = "drop"
                )
  return(list(Full=Data_Full,Summary=Data_Summary))
}

#' Plot clustering performance comparison
#'
#' @param df List. The output object from the \code{merge_simulation_results} function, 
#' containing \code{Full} and \code{Summary} data frames.
#' @param n Integer. The number of nodes.
#' @param d Numeric. The graph density.
#'
#' @return A \code{ggplot} object representing the performance grid across 
#' internal (\code{p_inside}) and external (\code{p_outside}) perturbation levels.
#' 
#' @import ggplot2
#' @export
plot_perf <- function(df, n, d) {
  # df is the output of merge_simulation_results function 
  p <- ggplot(df$Full,aes(x=Method, y=AMI,fill=Method)) +
    geom_boxplot(outlier.shape = NA)+scale_fill_manual(values=c("peachpuff","sandybrown","lightsalmon","lightsalmon3","tomato4","slategray2","steelblue1","steelblue3","steelblue4"))+
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

n<- 1000
d<- 0.005
df <- merge_simulation_results(path= "/Users/mchampion/Documents/GitHub/l1spectralclustering/",n,d)
g <- plot_perf(df, n,d)
g
