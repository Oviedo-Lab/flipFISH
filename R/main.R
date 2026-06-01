
# 1. Setup #############################################################################################################

# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2026

#' @useDynLib flipFISH, .registration = TRUE
#' @import Rcpp
#' @import RcppEigen 
#' @import ggplot2
NULL

# About: This script runs a simulation of "ground truth" RNA hybridization spots and noisy MERFISH barcode decoding. 
#   The aim is to compare the true precision (positive predictive value, PPV) to the confidence ratio (CR) of the barcodes 
#   while also assessing the use of blanks with count-based and CR-based thresholds as a FPQC metric. Simulations require 
#   a codebook, summary data statistics, and noise parameters. 

.onLoad <- function(libname, pkgname) {
#   # For plotting
#   title_size <- 18
#   axis_size <- 14
#   legend_size <- 14
#   dot_color <- "black"
#   filter_color <- "magenta"
#   filter_color_obs <- "magenta4"
#   filter_color2 <- "cyan"
#   filter_color2_obs <- "cyan4"
#   filter_style <- "solid"
#   retained_color <- "skyblue1"
#   retained_color_obs <- "skyblue4"
#   filtered_color <- "firebrick1"
#   filtered_color_obs <- "firebrick4"
#   blanks_color <- "gray"
#   blanks_color_obs <- "gray20"
#   
#   # For lazy parallelization 
#   cat("\nSetting up parallel processing")
#   if ( !is.element("package:parallel", search()) ) { # only try to load if not already loaded
#     if(!require(parallel)) {                         # if not installed, install
#       install.packages("parallel")                   # install if needed
#       library(parallel)                              # load
#     }
#   }
#   
#   # ... set core numbers
#   sys_name <- Sys.info()["sysname"]
#   core_num <- 1 # Default to 1 core for Windows, Windows doesn't support forking
#   if (sys_name == "Darwin" || sys_name == "Linux") core_num <- parallel::detectCores() - 2 
#   cat("\nNumber of cores to use:", core_num)
#   
#   # ... define helper function for running in parallel
#   runparallel <- function(input, funct, ...) {
#     results <- mclapply( 
#       X = input, FUN = funct, ...,
#       mc.preschedule = TRUE, mc.set.seed = TRUE,
#       mc.silent = FALSE, mc.cores = core_num,
#       mc.cleanup = TRUE, mc.allow.recursive = TRUE )
#     return(results)
#   }
#   
#   # Grab MERSCOPE codebook 
#   # ... provide your own codebook, or use the sample one
#   codebook <- read.csv("codebook_obs.csv", row.names = 1)
}

misreadQC <- function(
    STdata, # Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts", must have barcode names as row names
    codebook, # Codebook with row names as barcodes and columns as bits, must have barcode names as row names
    max_correctable_Hamming_distance,
  ) {
    
    # Get species names
    species_names <- rownames(STdata)
    # Align rows of codebook to STdata
    codebook <- codebook[species_names,]
    # Make blank mask
    blank_mask <- grepl("Blank", species_names)
    # Sort species by decreasing rates, with genes first and blanks second 
    gene_rate_order <- order(STdata$rates[!blank_mask], decreasing = TRUE)
    blank_rate_order <- order(STdata$rates[blank_mask], decreasing = TRUE)
    # Remake STdata and codebook with this order
    STdata <- rbind(
      STdata[!blank_mask,][gene_rate_order,],
      STdata[blank_mask,][blank_rate_order,]
    )
    codebook <- rbind(
      codebook[!blank_mask,][gene_rate_order,],
      codebook[blank_mask,][blank_rate_order,]
    )
    
    # Run misread QC algorithm with MCMCSA
    qc <- mQC(
      as.matrix(STdata), 
      as.matrix(codebook), 
      max_correctable_Hamming_distance,
      c(0.05, -(0.05 - 0.005)/100, 0.005), # step size, initial, slope, min
      c(0.1, -(0.1 - 0.01)/100, 0.01), # temp, initial, slope, min
      0.2,   # max_fr
      1e-6,   # ctol
      1000,    # n_steps
      5    # n_forks
    )
    
    # Make summary stats from qc results
    sum_names <- c("mean", "logmean", "ci_lower", "ci_upper")
    qc_names <- names(qc)
    bc_names <- qc_names[qc_names != "STdata" & qc_names != "fliprates"]
    bc_sum_names <- c()
    for (n in bc_names) {bc_sum_names <- c(bc_sum_names, paste0(n, "_", sum_names))}
    fr <- matrix(NA, nrow = ncol(qc$fliprates), ncol = 4)
    bc <- matrix(NA, nrow = ncol(qc$expected_corrected_counts), ncol = 4 * length(bc_names))
    colnames(fr) <- sum_names 
    colnames(bc) <- bc_sum_names
    rownames(bc) <- qc$STdata$species
    N_bits <- ncol(codebook)
    fr_names <- paste0("rate10_bit", seq_len(N_bits))
    fr_names <- c(fr_names, paste0("rate01_bit", seq_len(N_bits)))
    fr_names <- c(fr_names, paste0("corr_", seq_len(ncol(qc$fliprates) - 2*N_bits)))
    rownames(fr) <- fr_names
    for (p in qc_names) {
      if (p == "STdata") next
      p_means = colMeans(qc[[p]])
      p_log_means = colMeans(log(qc[[p]]))
      ci <- apply(qc[[p]], 2, quantile, probs = c(0.025, 0.975))
      if (p == "fliprates") {
        fr[,"mean"] <- p_means
        fr[,"logmean"] <- p_log_means
        fr[,"ci_lower"] <- ci[1,]
        fr[,"ci_upper"] <- ci[2,]
      } else {
        bc[,paste0(p, "mean")] <- p_means
        bc[,paste0(p, "logmean")] <- p_log_means
        bc[,paste0(p, "ci_lower")] <- ci[1,]
        bc[,paste0(p, "ci_upper")] <- ci[2,]
      }
    }
    qc[["fliprates_summary"]] <- fr
    qc[["bc_summary"]] <- bc
    
    return(qc)
    
  }

# Plot run QC (PPV cutoff)
plot_PPV_cutoff <- function(
    PPV,
    min_PPV = 0.8
  ) {
    
    # Sort by decreasing values
    PPV <- sort(PPV, decreasing = TRUE)
    PPV_cutoff <- rep("bad", length(PPV))
    PPV_cutoff[PPV > min_PPV] <- "good"
    df <- data.frame(x = seq_along(PPV), PPV = PPV, above_cutoff = PPV_cutoff)
    
    # Plot
    PPV_plot <- ggplot(df, aes(x = x, y = PPV)) +
      geom_point(aes(color = above_cutoff), size = 3) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
      scale_color_manual(values = c("good" = "blue", "bad" = "black")) +
      guides(color = "none") +
      labs(
        y = paste0("Estimated PPV"),
        x = paste0("Gene barcodes by PPV value"),
        title = "Estimated Precision (PPV) by Gene Barcode")
    # ^ Save as 1600 x 925
    
    return(PPV_plot)
    
  }

# Plot simulation QC
plot_counts_sorted <- function(
    counts_sorted,
    obs_counts_sorted,
    N_genes
  ) {
    
    bc_type_colors <- c(
      "Gene (sim)" = retained_color, "Blank (sim)" = blanks_color,
      "Gene (run)" = retained_color_obs, "Blank (run)" = blanks_color_obs
    )
    
    bc_type <- rep("Gene (sim)", length(counts_sorted))
    bc_type[seq_along(counts_sorted) > N_genes] <- "Blank (sim)"
    bc_type_obs <- rep("Gene (run)", length(obs_counts_sorted))
    bc_type_obs[seq_along(obs_counts_sorted) > N_genes] <- "Blank (run)"
    
    df_counts_sorted <- data.frame(
      index = c(seq_along(counts_sorted), seq_along(obs_counts_sorted)), 
      Count = c(counts_sorted, obs_counts_sorted),
      bc_type = c(bc_type, bc_type_obs),
      size = c(rep(3, length(counts_sorted)), rep(1.5, length(obs_counts_sorted)))
    )
    df_counts_sorted$Count[df_counts_sorted$Count == 0] <- 1
    df_counts_sorted$bc_type <- factor(
      df_counts_sorted$bc_type, 
      levels = c(
        "Gene (sim)", "Gene (run)", 
        "Blank (sim)", "Blank (run)"
      )
    )
    
    counts_sorted_plot <- ggplot(df_counts_sorted) +
      geom_point(size = df_counts_sorted$size, aes(x = index, y = Count, color = bc_type)) +
      scale_y_log10() +  # Log scale for better visibility
      scale_color_manual(values = bc_type_colors) +
      geom_vline(
        xintercept = mean(c(N_genes, N_genes+1)), 
        color = "black", 
        linetype = "dashed", 
        linewidth = 0.5
      ) +
      labs(title = "Simulated vs Run Counts", x = "Barcodes sorted by run count", y = "Spot count", color = "Count type") +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text.y = element_text(size = axis_size),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    return(counts_sorted_plot)
    
  }

# Plot simulation step information 
plot_MCMC_steps <- function(
    sim_summaries,
    burnin,
    FR_search
  ) {
    
    # Add step column
    sim_summaries$step <- c(1:nrow(sim_summaries))
    
    # Set v-line points
    burnin <- burnin + FR_search
    FR_search <- mean(c(FR_search, FR_search + 1))
    burnin <- mean(c(burnin, burnin + 1))
    vline_points <- c(FR_search, burnin)
   
    # plot sim_mse
    plot_mse <- ggplot(sim_summaries, aes(x = step, y = sim_mse)) +
      geom_line() +
      geom_vline(
        xintercept = vline_points, 
        color = "black", 
        linetype = "dashed", 
        linewidth = 0.5
      ) +
      labs(title = "MCMCSA step diagnostics", x = "step", y = "MSE (log)") +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text.y = element_text(size = axis_size),
        axis.text.x = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    # plot flip rates 
    df_flip <- rbind(
      data.frame(
        step = sim_summaries$step,
        FlipRate = sim_summaries$flip_rate_10_stip,
        Type = "1 to 0"),
      data.frame(
        step = sim_summaries$step,
        FlipRate = sim_summaries$flip_rate_01_stip,
        Type = "0 to 1")
    )
    
    # Plot
    plot_flip_rates <- ggplot(df_flip, aes(x = step, y = FlipRate, color = Type)) +
      geom_line() +
      geom_vline(
        xintercept = vline_points,
        color = "black",
        linetype = "dashed",
        linewidth = 0.5
      ) +
      labs(title = "MCMCSA Flip Rates", x = "step", y = "Flip Rate") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text.y = element_text(size = axis_size),
        axis.text.x = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    # plot bit noise correlations
    df_flip <- rbind(
      data.frame(
        step = sim_summaries$step,
        corr_param = sim_summaries$bit_noise_cor_sd_stip,
        Type = "sd"),
      data.frame(
        step = sim_summaries$step,
        corr_param = sim_summaries$bit_noise_cor_bias_stip,
        Type = "bias")
    )
    
    # Plot
    plot_corr_param <- ggplot(df_flip, aes(x = step, y = corr_param, color = Type)) +
      geom_line() +
      geom_vline(
        xintercept = vline_points,
        color = "black",
        linetype = "dashed",
        linewidth = 0.5
      ) +
      labs(title = "MCMCSA bit-noise correlation parameters", x = "step", y = "Value") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text.y = element_text(size = axis_size),
        axis.text.x = element_text(size = axis_size),
        legend.title = element_text(size = legend_size),
        legend.text = element_text(size = legend_size)
      )
    
    print(plot_mse)
    print(plot_flip_rates)
    print(plot_corr_param)
    
  }

