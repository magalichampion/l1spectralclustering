#' Merge simulation results across parameter grids
#'
#' @param path Character. The base directory path containing "data/" and "results/" subdirectories.
#' @param n Integer. The number of nodes in the graph.
#' @param k Numeric. Number of clusters.
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
merge_simulation_results <- function(path, n, k, p_inside_range = c(0.01,0.1,0.25,0.5), p_outside_range = c(0.01,0.1,0.25,0.5), method = c("l1Spectral", "Spectral", "regSpectral","robustSpectral", "VGAE", "ST_l1Spectral", "Spectrum", "Hybrid", "MCL")) {
  all_files <- list.files(paste0(path,"/results"), full.names = FALSE)
  
  data_full_list <- list()

  for (p_inside in p_inside_range){
    for (p_outside in p_outside_range){
      data_pattern <- paste0(path,"/data/data_n=", n,"_k=",sprintf("%g", k), "_p_inside=", sprintf("%.2f", p_inside), "_p_outside=", sprintf("%.2f", p_outside),".Rdata")
      if (!file.exists(data_pattern)) next
      load(data_pattern)
      
      for (current_method in method){
        file_pattern <- paste0("results_", current_method, 
                               "_n=", n, 
                               "_k=", sprintf("%g", k), 
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
              est_k <- NA
            } else {
              idx <- which(clus_est != 0) 
              est_k <- length(unique(clus_est[idx]))
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
              Estimated_k = est_k,  
              Method = current_method,
              p_inside = p_inside,
              p_outside = p_outside,
              n = n,
              k = k
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
                Mean_Estimated_k = mean(Estimated_k, na.rm = TRUE),
                .groups = "drop"
                )
  return(list(Full=Data_Full,Summary=Data_Summary))
}

#' Plot clustering performance comparison
#'
#' @param df List. The output object from the \code{merge_simulation_results} function, 
#' containing \code{Full} and \code{Summary} data frames.
#' @param n Integer. The number of nodes.
#' @param k Numeric. Number of clusters.
#'
#' @return A \code{ggplot} object representing the performance grid across 
#' internal (\code{p_inside}) and external (\code{p_outside}) perturbation levels.
#' 
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
plot_perf <- function(df, n, k) {

  # 1. Reshape data for l1 methods
  l1_methods <- c("l1Spectral", "ST_l1Spectral")
  df_plot <- df$Full %>% mutate(Metric = "Standard")
  df_corrected <- df_plot %>%
    filter(Method %in% l1_methods) %>%
    mutate(AMI = AMI_corrected, Metric = "Corrected")
  
  df_final <- bind_rows(df_plot, df_corrected) %>%
    mutate(Metric = factor(Metric, levels = c("Standard", "Corrected")))
  
  # 2. Colors & categories
  methods_non  <- unique(df$Full$Method[df$Full$Category == "Non self-tuned"])
  methods_self <- unique(df$Full$Method[df$Full$Category == "Self-tuned"])
  
  warm_refined <- c("#FAD7A0", "#F8C471", "#E59866", "#D35400", "#A04000")
  cool_refined <- c("#AED6F1", "#5DADE2", "#2E86C1", "#1B4F72")
  final_palette <- c(warm_refined[1:length(methods_non)], 
                     cool_refined[1:length(methods_self)])
  
  # 3. Main plot
  p <- ggplot(df_final, aes(x = Method, y = AMI, fill = Method)) +
    geom_boxplot(
      aes(alpha = Metric), 
      outlier.shape = NA, 
      linewidth = 0.4, 
      fatten = 1.5,
      position = position_dodge(width = 0.8)
    ) +
    scale_alpha_manual(values = c("Standard" = 1, "Corrected" = 0.35)) +
    scale_fill_manual(values = final_palette) +
    facet_grid(p_outside ~ p_inside, labeller = labeller(
      p_inside = function(x) paste0("p_in = ", x),
      p_outside = function(x) paste0("p_out = ", x)
    )) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  p <- p + guides(
    fill = guide_legend(
      title = "Method Category (Row 1: Non Self-Tuned | Row 2: Self-Tuned)",
      nrow = 2, byrow = TRUE, order = 1
    ),
    alpha = guide_legend(
      title = "AMI Type:",
      # Using gray50 makes the transparency visible without picking a specific method's color
      override.aes = list(fill = "gray50"), 
      order = 2
    )
  )
  
  p <- p + geom_vline(xintercept = length(methods_non) + 0.5, linetype = "dotted")
  
  return(p + labs(title = paste0("n = ", n, " and k = ", k), y = "AMI"))
}

#' Plot estimated number of clusters
#'
#' @param df List. The output object from the \code{merge_simulation_results} function, 
#' containing \code{Full} and \code{Summary} data frames.
#' @param n Integer. The number of nodes.
#' @param k Numeric. True number of clusters.
#'
#' @return A \code{ggplot} object representing the estimated number of clusters across 
#' internal (\code{p_inside}) and external (\code{p_outside}) perturbation levels.
#' 
#' @import ggplot2
#' @export
plot_n_cluster <- function(df, n, k){
  df_plot <- df$Full %>% filter(Method %in% c("ST_l1Spectral", "Spectrum", "Hybrid", "MCL")) %>%
    mutate(Method = factor(Method))
  
  # 1. Define order and palette
  methods_order <- c("ST_l1Spectral", "Spectrum", "Hybrid", "MCL")
  cool_refined  <- c("#AED6F1", "#5DADE2", "#2E86C1", "#1B4F72")
  names(cool_refined) <- methods_order
  
  # 2. Main plot
  g <- ggplot(df_plot, aes(x = Method, y = Estimated_k, fill = Method)) +
    
    # Layer 1: Black jittered dots
    geom_jitter(
      color = "black",      
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.4,          
      size = 0.4
    ) + 
    
    # Layer 2: Boxplot with colored fill and black borders/medians 
    geom_boxplot(
      color = "black",      # Black borders and median line
      outlier.shape = NA, 
      linewidth = 0.5, 
      fatten = 2,           
      position = position_dodge(width = 0.8)
    ) +
    
    # Apply the blue palette to the fill of the boxes
    scale_fill_manual(values = cool_refined) +
    
    # Facet layout 
    facet_grid(p_outside ~ p_inside, labeller = labeller(
      p_inside = function(x) paste0("p_in = ", x),
      p_outside = function(x) paste0("p_out = ", x)
    )) +
    
    # Reference line for ground truth 
    geom_hline(yintercept = k, linetype = "dashed", color = "black", alpha = 0.5) +
    
    # Styling and theme
    scale_y_continuous(breaks = c(1, 10, 20, 30)) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = expression(Estimated~Number~of~Clusters),
      x = "Method"
    )
  
  # 4. Legend styling
  g <- g + guides(
    fill = guide_legend(nrow = 1)
  )
}