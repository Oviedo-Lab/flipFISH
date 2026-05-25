
# 1. Setup #############################################################################################################

# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2026

#' @useDynLib FISHmQC, .registration = TRUE
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




# 5. Functions to run and analyze simulation ###########################################################################

# Function to run individual simulation
run_simulation <- function(
    codebook,
    summary_stats_genes,
    summary_stats_blanks,
    flip_rate_10,                           # Probability of a 0 being flipped to a 1, either scalar or vector of length equal to number of bits
    flip_rate_01,                           # Probability of a 1 being flipped to a 0, either scalar or vector of length equal to number of bits
    max_correctable_Hamming_distance,       # For the automatic error correction, assumes an "extended" Hamming code with min 4 bits between all barcodes
    bit_noise_correlations_bias = 0,        # For correlating noise in the bit reads
    bit_noise_correlations_sd = 1,          # For correlating noise in the bit reads
    bit_lum_noise_correlations = NULL,      # For stipulating correlations between bits (bit bias and sd will be ignored)
    N_cells = 1e3,                          # Number of cells to simulate
    blank_weight = 1,
    maintain_gene_identity = TRUE, ## actual sampling (and burn-in) done with this set to "TRUE"; only the initial fr search steps use "FALSE". 
    return_flip_rates = FALSE,
    verbose = FALSE
  ) {
   
    # Extract summary stats
    observed_gene_rates <- summary_stats_genes$rates
    observed_gene_variance <- summary_stats_genes$variance
    observed_blank_rates <- summary_stats_blanks$rates
    observed_blank_variance <- summary_stats_blanks$variance
    
    # Extract simulation parameters
    N_bits <- ncol(codebook)
    N_barcodes <- nrow(codebook)
    N_blanks <- length(observed_blank_rates)
    N_genes <- length(observed_gene_rates)
    mean_Hamming_weight <- mean(rowSums(codebook))
    # ... calculate number of bit-pair correlations
    N_corrs <- 0
    for (bits_left in N_bits:1) N_corrs <- N_corrs + bits_left
    N_corrs <- N_corrs - N_bits
   
    # Report simulation parameters
    if (verbose) {
      cat("\nSimulation parameters set:")
      cat("\n\tNumber of cells:", N_cells)
      cat("\n\tNumber of bits:", N_bits)
      cat("\n\tNumber of barcodes:", N_barcodes)
      cat("\n\tNumber of genes:", N_genes)
      cat("\n\tNumber of blanks:", N_blanks)
      cat("\n\tNumber of bit-pair correlations:", N_corrs)
      cat("\n\tMean Hamming weight:", mean_Hamming_weight)
      if (is.null(bit_lum_noise_correlations)) {
        cat("\n\tBit noise correlations bias:", bit_noise_correlations_bias)
        cat("\n\tBit noise correlations sd:", bit_noise_correlations_sd)
      } else {
        cat("\n\tBit noise correlations bias:", mean(bit_lum_noise_correlations))
        cat("\n\tBit noise correlations sd:", sd(bit_lum_noise_correlations))
      }
      cat("\n\tFlip rate 1 -> 0:", mean(flip_rate_10))
      cat("\n\tFlip rate 0 -> 1:", mean(flip_rate_01))
      cat("\n\tMax correctable Hamming distance:", max_correctable_Hamming_distance)
    }
    
    # Make generation matrix for codebook, representing the expected luminance values for each bit in each barcode
    # ... regenerated with random (but properly valenced) numbers each simulation
    if (verbose) cat("\nMaking codebook generation matrix")
    codebook_generation_matrix <- make_codebook_generation_matrix(codebook)
    
    # Make blank mask from row names
    if (verbose) cat("\nMaking blank mask from codebook row names")
    blank_mask <- rep(FALSE, N_barcodes)
    if (length(rownames(codebook) == N_barcodes)) blank_mask <- grepl("Blank", rownames(codebook))
    if (!any(blank_mask)) {
      cat("\nNo blanks detected in codebook row names, setting blanks randomly. This could be a problem.")
      blank_mask[sample(1:nrow(codebook), size = N_blanks, replace = FALSE)] <- TRUE
    }
    
    # Set the correlation matrix for luminance noise
    # ... correlations set randomly each simulation
    if (verbose) cat("\nSetting correlation matrix for luminance noise")
    Sigma <- set_luminance_noise_correlation_matrix(
      noise_correlation_bias = bit_noise_correlations_bias, 
      noise_correlation_sd = bit_noise_correlations_sd,
      lum_noise_correlations = bit_lum_noise_correlations,
      N_bits = N_bits,
      N_corrs = N_corrs
    )
    
    # Predict decoding rate
    expected_flip_rate <- (mean(flip_rate_10)*mean_Hamming_weight + mean(flip_rate_01)*(N_bits - mean_Hamming_weight)) / N_bits
    expected_decoding_rate <- ppois(q = max_correctable_Hamming_distance/2, lambda = expected_flip_rate * N_bits)
    if (verbose) cat("\nPredicted decoding rate:", expected_decoding_rate)
    
    # Set true counts
    # ... set randomly each simulation, based around the summary statistics of the observed data
    if (verbose) cat("\nSetting true counts for each gene and blank")
    true_spot_info <- set_true_spot_info(
      observed_gene_rates/expected_decoding_rate, 
      observed_gene_variance, 
      blank_mask,
      N_cells,
      N_barcodes,
      N_blanks,
      N_genes
    )
    
    # Set noise based on requested flip rates
    # ... set randomly each simulation, with expected flip rates
    if (verbose) cat("\nSetting simulated noise (sd) needed in each bit for stipulated flip rates")
    codebook_noise <- set_noise_by_flip_rates(
      codebook_generation_matrix = codebook_generation_matrix, 
      flip_rate_10 = flip_rate_10, 
      flip_rate_01 = flip_rate_01
    )
    
    # Generate simulated spots
    if (verbose) cat("\nGenerating simulated spots")
    spots <- runparallel(
      1:N_cells, 
      generate_simulated_spots, 
      noise = codebook_noise,
      bc_counts = true_spot_info$barcode_counts_true,
      codebook_gen_mat = codebook_generation_matrix,
      corr_mat = Sigma,
      N_bits = N_bits,
      decode_and_label = !verbose && !return_flip_rates, 
      codebook = codebook, 
      max_correctable_Hamming_distance = max_correctable_Hamming_distance
    )
    
    # ... initialize flip rates
    # ... computationally, very expensive to check, only check if verbose
    flip_rates <- c(NA, NA)
    names(flip_rates) <- c("flip_rate_10", "flip_rate_01")
    
    if (verbose) {
      
      # Check luminance noise correlations
      cat("\nChecking luminance noise correlations by bit")
      # ... find noise by bit by cell
      luminance_noise <- runparallel(
        1:N_cells, 
        find_noise,
        spots = spots,
        bc_counts = true_spot_info$barcode_counts_true,
        codebook_gen_mat = codebook_generation_matrix,
        N_bits = N_bits
      )
      # ... collapse into single matrix of noise values (rows as spots, columns as bits)
      luminance_noise <- do.call(rbind, luminance_noise)
      # ... calculate correlations in noise between bits
      corrs_all <- matrix(NA, nrow = N_corrs, ncol = 2)
      idx <- 0
      for (bi in 1:(N_bits-1)) {
        for (bj in (bi+1):N_bits) {
          idx <- idx + 1
          noise_correlation <- cor(
            luminance_noise[,bi], 
            luminance_noise[,bj] 
          )
          expected_correlation <- Sigma[bi,bj]
          corrs_all[idx,] <- c(noise_correlation, expected_correlation)
        }
      }
      # ... print results
      cat("\nCorrelation between observed and expected luminance-noise correlations:", 
          cor(corrs_all[,1], corrs_all[,2], use = "pairwise.complete.obs"))
      
      # Decode spots 
      cat("\nDecoding spots")
      spots <- runparallel(
        1:N_cells, 
        decode_spots, 
        spots = spots,
        threshold = 0
      )
      
      # Check unique correctable error sizes
      cat("\nChecking unique correctable error sizes")
      corrected_distances <- check_correctable_distances(
        spots = spots, 
        codebook = codebook,
        max_correctable_Hamming_distance = max_correctable_Hamming_distance
      )
      cat("\nUnique correctable error sizes:", 
          paste(corrected_distances, collapse = ", "))
      
      # Find the ground truth behind the simulated spots
      cat("\nGenerating simulated mRNA molecules responsible for those spots (ground truth)")
      RNA <- runparallel(
        1:N_cells, 
        generate_RNA,
        bc_counts = true_spot_info$barcode_counts_true,
        codebook = codebook,
        N_bits = N_bits
      )
      
      # Check flip rates
      cat("\nChecking bit-flip rates")
      flip_rates <- check_flip_rates(
        RNA, 
        spots, 
        verbose = verbose
      )
      
      # Find uncorrected bc and apply read corrections
      if (verbose) cat("\nFinding uncorrected barcode labels and applying read corrections")
      spots <- runparallel(
        1:N_cells, 
        apply_spot_correction, 
        spots = spots, 
        codebook = codebook, 
        max_correctable_Hamming_distance = max_correctable_Hamming_distance
      )
      
    } else if (return_flip_rates) {
      
      # Decode spots 
      spots <- runparallel(
        1:N_cells, 
        decode_spots, 
        spots = spots,
        threshold = 0
      )
      
      # Find the ground truth behind the simulated spots
      RNA <- runparallel(
        1:N_cells, 
        generate_RNA,
        bc_counts = true_spot_info$barcode_counts_true,
        codebook = codebook,
        N_bits = N_bits
      )
      
      # Check flip rates
      flip_rates <- check_flip_rates(
        RNA, 
        spots, 
        verbose = verbose
      )
      
      # Find uncorrected bc and apply read corrections
      spots <- runparallel(
        1:N_cells, 
        apply_spot_correction, 
        spots = spots, 
        codebook = codebook, 
        max_correctable_Hamming_distance = max_correctable_Hamming_distance
      )
      
    }
    
    spots_corrected_bc_labels <- lapply(spots, `[[`, "spot_labels")
    spots_bc_labels <- lapply(spots, `[[`, "spot_labels_uncorrected")
    spot_decoding_rate <- find_decoded_spot_rate(spots_bc_labels)
    spot_decoding_rate_corrected <- find_decoded_spot_rate(spots_corrected_bc_labels)
    if (verbose) {
      cat("\nSpot decoding rate, no error correction:", spot_decoding_rate)
      cat("\nSpot decoding rate, with error correction:", spot_decoding_rate_corrected)
    }
    
    # Compute PPV
    if (verbose) cat("\nComputing PPV")
    PPV_corrected_decoding <- compute_PPV_fast(
      true_spot_info$barcodes_by_cell_true, 
      spots_corrected_bc_labels,
      N_barcodes
    )
    PPV_corrected_decoding <- PPV_corrected_decoding[!blank_mask]
    
    # Extract barcode counts from spot labels 
    if (verbose) cat("\nExtracting barcode counts from spot labels")
    barcode_counts_est <- find_barcode_counts(spots_bc_labels, N_barcodes)
    barcode_counts_est_corrected <- find_barcode_counts(spots_corrected_bc_labels, N_barcodes)
    # ... extract results 
    total_barcode_counts_est <- colSums(barcode_counts_est, na.rm = TRUE)
    total_barcode_counts_est_genes <- total_barcode_counts_est[!blank_mask] 
    total_barcode_counts_est_corrected <- colSums(barcode_counts_est_corrected, na.rm = TRUE)
    total_barcode_counts_est_corrected_genes <- total_barcode_counts_est_corrected[!blank_mask] 
    total_barcode_counts_est_corrected_blanks <- total_barcode_counts_est_corrected[blank_mask] 
    max_blank_count <- max(total_barcode_counts_est_corrected_blanks)
    sim_gene_order <- order(total_barcode_counts_est_corrected_genes, decreasing = TRUE)
    sim_blank_order <- order(total_barcode_counts_est_corrected_blanks, decreasing = TRUE)
    
    # ... while we're at it, make sorted vector for the real data
    obs_counts_genes <- observed_gene_rates * N_cells
    obs_counts_blanks <- observed_blank_rates * N_cells
    obs_counts <- rep(NA, N_barcodes)
    obs_counts[!blank_mask] <- obs_counts_genes
    obs_counts[blank_mask] <- obs_counts_blanks
    max_blank_count_obs <- max(obs_counts_blanks)
    obs_gene_order <- order(obs_counts_genes, decreasing = TRUE)
    obs_blank_order <- order(obs_counts_blanks, decreasing = TRUE)
    obs_counts_sorted <- c(
      obs_counts_genes[obs_gene_order], 
      obs_counts_blanks[obs_blank_order] # use same order as in simulation so we can compare relative bit-flip rates
    )
    
    if (maintain_gene_identity) {
      counts_sorted <- c(
        total_barcode_counts_est_corrected_genes[obs_gene_order], 
        total_barcode_counts_est_corrected_blanks[obs_blank_order]
      )
    } else {
      counts_sorted <- c(
        total_barcode_counts_est_corrected_genes[sim_gene_order], 
        total_barcode_counts_est_corrected_blanks[sim_blank_order] 
      )
    }
    
    # Find mse of the log of the data
    gene_weight <- 1 + (N_blanks/N_genes) * (1 - blank_weight)
    weights <- c(rep(gene_weight, N_genes), rep(blank_weight, N_blanks))
    sim_se <- (log(counts_sorted + 1) - log(obs_counts_sorted + 1))^2
    sim_mse <- mean(sim_se * weights)
    if (verbose) cat("\nWeighted mean squared log-error:", sim_mse)
    
    # Compute confidence ratio for barcode counts
    if (verbose) cat("\nComputing confidence ratio")
    CR <- total_barcode_counts_est / total_barcode_counts_est_corrected
    CR[total_barcode_counts_est_corrected == 0] <- 0
    if (any(is.na(CR) | is.infinite(CR) | is.nan(CR) | isTRUE(CR <= 0))) {
      stop("Confidence ratio contains invalid values (NA, Inf, NaN, or <= 0).")
    }
    # ... extract results
    CR_genes <- CR[!blank_mask]
    CR_blanks <- CR[blank_mask]
    max_blank_CR <- max(CR_blanks)
    CR_sorted <- c(
      CR_genes[order(CR_genes, decreasing = TRUE)], 
      CR_blanks[order(CR_blanks, decreasing = TRUE)]
    )
    
    # Find relationship between CR and PPV 
    # ... fit linear model with intercept = 0
    CR_PPV_lm <- lm(PPV_corrected_decoding ~ CR_genes + 0)  
    # ... extract slope
    CR_PPV_lm_slope <- coef(CR_PPV_lm)[["CR_genes"]]
    # ... compute predictions and MSE
    predictions <- CR_PPV_lm_slope * CR_genes
    mse <- mean((PPV_corrected_decoding - predictions)^2)
    if (verbose) {
      cat("\nLinear model of PPV ~ CR + 0:")
      cat("\n\tSlope:", CR_PPV_lm_slope)
      cat("\n\tMSE:", mse)
    }
    
    # Extract cutoff stats
    # ... assume that PPV is normally distributed for genes below the max-count blank cutoff 
    retained_gene_mask <- total_barcode_counts_est_corrected_genes > max_blank_count
    PPV_mean_filtered <- mean(PPV_corrected_decoding[!retained_gene_mask])
    PPV_sd_filtered <- sd(PPV_corrected_decoding[!retained_gene_mask])
    # ... assume that PPV is exponentially distributed for genes above the max-count blank cutoff
    PPV_mean_retained <- mean(PPV_corrected_decoding[retained_gene_mask])
    # ... get percent of genes filtered by CR cutoff
    ratio_genes_CR_filtered <- sum(CR_genes <= max_blank_CR) / length(CR_genes)
    
    # Gather results
    bit_noise_correlations_bias_out <- bit_noise_correlations_bias
    bit_noise_correlations_sd_out <- bit_noise_correlations_sd
    if (!is.null(bit_lum_noise_correlations)) {
      bit_noise_correlations_bias_out <- mean(bit_lum_noise_correlations)
      bit_noise_correlations_sd_out <- sd(bit_lum_noise_correlations)
    }
    sim_summary <- c(
      CR_PPV_lm_slope, 
      mse,
      sim_mse,
      PPV_mean_filtered,
      PPV_sd_filtered,
      PPV_mean_retained,
      ratio_genes_CR_filtered,
      mean(flip_rate_10),
      flip_rates[1], 
      mean(flip_rate_01),
      flip_rates[2],
      bit_noise_correlations_bias_out, bit_noise_correlations_sd_out,
      spot_decoding_rate, spot_decoding_rate_corrected
      )
    names(sim_summary) <- c(
      "CR_PPV_slope", 
      "CR_PPV_mse",
      "sim_mse",
      "PPV_mean_filtered", 
      "PPV_sd_filtered",
      "PPV_mean_retained",
      "ratio_genes_CR_filtered",
      "flip_rate_10_stip",
      names(flip_rates)[1], 
      "flip_rate_01_stip",
      names(flip_rates)[2],
      "bit_noise_cor_bias_stip", "bit_noise_cor_sd_stip",
      "spot_decoding_rate", "spot_decoding_rate_corrected"
      )
    
    # Each row is a simulated gene: 
    gene_summary <- data.frame(
      Count = total_barcode_counts_est_genes,
      Count_corrected = total_barcode_counts_est_corrected_genes,
      CR = CR_genes,
      PPV = PPV_corrected_decoding
    )
    rownames(gene_summary) <- rownames(codebook)[!blank_mask]
    
    results <- list(
      gene_summary = gene_summary,
      counts_sorted = data.frame(sim = counts_sorted, obs = obs_counts_sorted),
      sim_summary = sim_summary
    )
    
    if (verbose) {
      
      blank_line_width <- 1.0
      observed_gene_counts <- summary_stats_genes$counts
      observed_blank_counts <- summary_stats_blanks$counts
      
      # barcode count distributions
      df_obs_est <- data.frame(x = log(c(observed_gene_counts, observed_blank_counts)+1))
      df_sim_est <- data.frame(x = log(colSums(barcode_counts_est_corrected)+1))
      df_sim_true <- data.frame(x = log(colSums(true_spot_info$barcode_counts_true)+1))
      
      cat("\nMean spots per cell:",
          "\n\tObserved, decoded:", mean(c(observed_gene_rates, observed_blank_rates)), 
          "\n\tSimulated, decoded:", mean(colSums(barcode_counts_est_corrected)/N_cells), 
          "\n\tSimulated, true:", mean(colSums(true_spot_info$barcode_counts_true)/N_cells))
      
      plot_obs_est <- ggplot(df_obs_est, aes(x)) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black") +
        labs(title = "Real decoded counts", x = "Barcode count", y = NULL) + 
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = title_size*0.8),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      plot_sim_est <- ggplot(df_sim_est, aes(x)) +
        geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
        labs(title = "Simulated decoded counts", x = "Barcode count", y = NULL) + 
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = title_size*0.8),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      plot_sim_true <- ggplot(df_sim_true, aes(x)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        labs(title = "Simulated ground-truth counts", x = "Barcode count", y = NULL) + 
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = title_size*0.8),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      grid.arrange(
        plot_obs_est, plot_sim_est, plot_sim_true, 
        nrow = 1, 
        top = textGrob(
          "Barcode Count Frequency", 
          gp = gpar(fontsize = title_size)
        ) 
      )
      
      bc_type_colors <- c(
        "Retained" = retained_color, "Filtered" = filtered_color, "Blanks" = blanks_color,
        "Retained_obs" = retained_color_obs, "Filtered_obs" = filtered_color_obs, "Blanks_obs" = blanks_color_obs
      )
      
      # CR filter
      bc_type <- rep("Retained", length(CR_sorted))
      bc_type[CR_sorted <= max_blank_CR] <- "Filtered"
      bc_type[seq_along(CR_sorted) > length(CR_genes)] <- "Blanks"
      df_CR_sorted <- data.frame(
        index = seq_along(CR_sorted), 
        CR = CR_sorted,
        bc_type = bc_type
      )
      df_CR_sorted$CR[df_CR_sorted$CR == 0] <- min(df_CR_sorted$CR[df_CR_sorted$CR > 0]) * 0.01
      df_CR_sorted$bc_type <- factor(df_CR_sorted$bc_type, levels = c("Retained", "Filtered", "Blanks"))
      plot_CR_sorted <- ggplot(df_CR_sorted, aes(x = index, y = CR, color = bc_type)) +
        geom_point(size = 3) +
        scale_y_log10() +  # Log scale for better visibility
        scale_color_manual(values = bc_type_colors) +
        geom_vline(
          xintercept = mean(c(length(CR_genes), length(CR_genes)+1)), 
          color = "black", 
          linetype = "dashed", 
          linewidth = 0.5
        ) +
        geom_vline(
          xintercept = mean(c(min(which(CR_sorted <= max_blank_CR)), min(which(CR_sorted <= max_blank_CR)) + 1)), 
          color = filter_color, 
          linetype = filter_style, 
          linewidth = blank_line_width
        ) +
        labs(title = "Max-CR filter", x = "Barcodes sorted by CR", y = "Confidence ratio (CR)", color = "Barcode status") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text.y = element_text(size = axis_size),
          axis.text.x = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      # Count filter
      bc_type <- rep("Retained", length(counts_sorted))
      bc_type[counts_sorted <= max_blank_count] <- "Filtered"
      bc_type[seq_along(counts_sorted) > length(total_barcode_counts_est_corrected_genes)] <- "Blanks"
      bc_type_obs <- rep("Retained_obs", length(obs_counts_sorted))
      bc_type_obs[obs_counts_sorted <= max_blank_count_obs] <- "Filtered_obs"
      bc_type_obs[seq_along(obs_counts_sorted) > length(obs_counts_genes)] <- "Blanks_obs"
      df_counts_sorted <- data.frame(
        index = c(seq_along(counts_sorted), seq_along(obs_counts_sorted)), 
        Count = c(counts_sorted, obs_counts_sorted),
        bc_type = c(bc_type, bc_type_obs),
        data_type = c(rep("Simulated", length(counts_sorted)), rep("Observed", length(obs_counts_sorted))),
        size = c(rep(3, length(counts_sorted)), rep(1.5, length(obs_counts_sorted)))
      )
      df_counts_sorted$Count[df_counts_sorted$Count == 0] <- 1
      df_counts_sorted$bc_type <- factor(
        df_counts_sorted$bc_type, 
        levels = c(
          "Retained", "Retained_obs", "Filtered",
          "Filtered_obs", "Blanks", "Blanks_obs"
          )
        )
      plot_counts_sorted <- ggplot(df_counts_sorted) +
        geom_point(size = df_counts_sorted$size, aes(x = index, y = Count, color = bc_type)) +
        scale_y_log10() +  # Log scale for better visibility
        scale_color_manual(values = bc_type_colors) +
        geom_vline(
          xintercept = mean(c(min(which(obs_counts_sorted <= max_blank_count_obs)), min(which(obs_counts_sorted <= max_blank_count_obs)) + 1)), 
          color = filter_color_obs, 
          linetype = filter_style, 
          linewidth = blank_line_width/2
        ) +
        geom_vline(
          xintercept = mean(c(length(total_barcode_counts_est_corrected_genes), length(total_barcode_counts_est_corrected_genes)+1)), 
          color = "black", 
          linetype = "dashed", 
          linewidth = 0.5
        ) +
        geom_vline(
          xintercept = mean(c(min(which(counts_sorted <= max_blank_count)), min(which(counts_sorted <= max_blank_count)) + 1)), 
          color = filter_color, 
          linetype = filter_style, 
          linewidth = blank_line_width
        ) +
        labs(title = "Max-count filter", x = "Barcodes sorted by count", y = "Total counts per barcode", color = "Barcode status") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text.y = element_text(size = axis_size),
          axis.text.x = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      print(plot_CR_sorted)
      print(plot_counts_sorted)
      
      # # Extract the legend from one plot
      # g_legend <- function(a.gplot) {
      #   tmp <- ggplotGrob(a.gplot)
      #   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      #   legend <- tmp$grobs[[leg]]
      #   return(legend)
      # }
      # 
      # legend <- g_legend(plot_counts_sorted)
      # 
      # plot_CR_sorted <- plot_CR_sorted + theme(legend.position = "none")
      # plot_counts_sorted <- plot_counts_sorted + theme(legend.position = "none")
      # 
      # grid.arrange(
      #   arrangeGrob(plot_counts_sorted, plot_CR_sorted, ncol = 2),
      #   legend,
      #   ncol = 1,
      #   heights = c(9, 1)  # Adjust as needed
      # )
      
      # CR vs PPV (genes only)
      ret <- rep("no", length(CR_genes))
      ret[CR_genes > max_blank_CR] <- "yes"
      df_CR_filter <- data.frame(
        CR = CR_genes,
        PPV = PPV_corrected_decoding,
        Retained = ret
      )
      plot_CR_filter <- ggplot(df_CR_filter, aes(x = CR, y = PPV, color = Retained)) +
        geom_point() +
        scale_color_manual(values = c("no" = "gray", "yes" = dot_color)) +
        geom_vline(
          xintercept = max_blank_CR, 
          color = filter_color, 
          linetype = filter_style, 
          linewidth = blank_line_width
        ) +
        ylim(c(0, 1)) +  # Limit y-axis to [0, 1]
        xlim(c(0, 1)) +
        labs(title = "Max blank CR", x = "Confidence ratio (CR)", y = NULL) + 
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = title_size*0.8),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      # Total Barcode Counts vs PPV (genes only)
      ret <- rep("no", length(total_barcode_counts_est_corrected_genes))
      ret[total_barcode_counts_est_corrected_genes > max_blank_count] <- "yes"
      df_count_filter <- data.frame(
        Count = total_barcode_counts_est_corrected_genes,
        PPV = PPV_corrected_decoding,
        Retained = ret
      )
      plot_count_filter <- ggplot(df_count_filter, aes(x = Count, y = PPV, color = Retained)) +
        geom_point() +
        scale_x_log10() +  # Log scale for better visibility
        scale_color_manual(values = c("no" = "gray", "yes" = dot_color)) +
        geom_vline(
          xintercept = max_blank_count, 
          color = filter_color, 
          linetype = filter_style, 
          linewidth = blank_line_width
        ) +
        labs(title = "Max blank count", x = "Total counts per barcode", y = NULL) + 
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = title_size*0.8),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size)
        )
      
      filter_comp <- arrangeGrob(
        plot_count_filter, plot_CR_filter,  
        nrow = 1, 
        top = textGrob(
          "Comparison of blank-based filters for genes", 
          gp = gpar(fontsize = title_size)
          ) 
        )
      
      grid.arrange(
        textGrob("PPV", rot = 90, gp = gpar(fontsize = axis_size)),
        filter_comp,
        ncol = 2,
        widths = unit.c(unit(1, "lines"), unit(1, "null"))
      )
      
    }
    
    return(results)
    
  }

# Estimate FPQC metrics for data using Markov chain Monte Carlo simulated annealing
run_MCMCSA <- function(
    codebook,
    mouse_id,                      # m1, m2, m3, etc...
    resamples,                     # Number of steps to take in random walk
    burnin = NULL,
    FR_search = NULL,
    max_step_size = 0.2,         # Size of each step in the random walk
    min_step_size_divider = 10,   # Scale step-size down linearly to max_step_size/min_step_size_divider over burnin
    max_temperature = 0.2,     # Temperature for the MCMC simulation
    min_temperature_scale = 0.5, # Scale temperature down linearly to max_temperature/(min_step_size_divider/min_temperature_scale) over burnin
    max_blank_weight = NULL,
    min_blank_weight = NULL,
    n_tracker_updates = 100
  ) {

    cat("\n\nEstimating FPQC metrics with Markov chain Monte Carlo simulated annealing")

    # Get summary stats
    summary_stats_genes <- summary_stats(paste0(mouse_id, "_genes"))
    summary_stats_blanks <- summary_stats(paste0(mouse_id, "_blanks"))

    # Initialize step counter
    step <- 0
    # Initialize call counter
    ctr <- 0

    # Extract simulation parameters
    N_bits <- ncol(codebook)
    N_blanks <- length(summary_stats_blanks$rates)
    N_genes <- length(summary_stats_genes$rates)
    # ... calculate number of bit-pair correlations
    N_corrs <- 0
    for (bits_left in N_bits:1) N_corrs <- N_corrs + bits_left
    N_corrs <- N_corrs - N_bits

    # Check max_blank_weight
    max_blank_weight_limit <- N_genes/N_blanks + 1
    if (is.null(max_blank_weight)) {
      max_blank_weight <- max_blank_weight_limit
    } else if (max_blank_weight > max_blank_weight_limit) {
      max_blank_weight <- max_blank_weight_limit
      cat("\n\nRequested max blank weight too high, resetting to limit:", max_blank_weight_limit)
    }

    # Set up steps
    if (is.null(FR_search)) FR_search <- as.integer(resamples/2)
    if (is.null(burnin)) burnin <- as.integer(resamples/2)
    n_steps <- FR_search + burnin + resamples

    # Set up schedules
    # ... set mins
    min_step_size <- max_step_size/min_step_size_divider
    min_temperateure <- max_temperature/(min_step_size_divider/min_temperature_scale)
    if (is.null(min_blank_weight)) {
      min_blank_weight <- 2
    } else if (min_blank_weight < 1) {
      min_blank_weight <- 1
    }

    # ... set flip-rate search schedules
    step_size_schedule <- seq(from = max_step_size, to = min_step_size, length.out = FR_search)
    temperature_schedule <- seq(from = max_temperature, to = min_temperateure, length.out = FR_search)
    blank_weight_schedule <- seq(from = max_blank_weight, to = min_blank_weight, length.out = FR_search)
    n_retained_bits <- as.integer(seq(from = N_bits - 1, to = 0, length.out = FR_search))
    n_retained_corrs <- as.integer(seq(from = N_corrs - 1, to = 0, length.out = FR_search))

    # ... set burnin
    step_size_schedule <- c(step_size_schedule, seq(from = max_step_size, to = min_step_size, length.out = burnin))
    temperature_schedule <- c(temperature_schedule, seq(from = max_temperature, to = min_temperateure, length.out = burnin))
    blank_weight_schedule <- c(blank_weight_schedule, seq(from = max_blank_weight, to = min_blank_weight, length.out = burnin))
    n_retained_bits <- c(n_retained_bits, as.integer(seq(from = N_bits - 1, to = 0, length.out = burnin)))
    n_retained_corrs <- c(n_retained_corrs, as.integer(seq(from = N_corrs - 1, to = 0, length.out = burnin)))

    # ... set sampling
    step_size_schedule <- c(step_size_schedule, rep(min_step_size, resamples))
    temperature_schedule <- c(temperature_schedule, rep(min_temperateure, resamples))
    blank_weight_schedule <- c(blank_weight_schedule, rep(min_blank_weight, resamples))
    n_retained_bits <- c(n_retained_bits, rep(0, resamples))
    n_retained_corrs <- c(n_retained_corrs, rep(0, resamples))

    # Initialize structures to hold results
    sim_summaries <- as.data.frame(matrix(NA, nrow = n_steps, ncol = 15))
    PPV_genes <- as.data.frame(matrix(NA, nrow = n_steps, ncol = N_genes))
    Counts_sorted_sim <- as.data.frame(matrix(NA, nrow = n_steps, ncol = N_genes + N_blanks))
    Counts_sorted_obs <- as.data.frame(matrix(NA, nrow = n_steps, ncol = N_genes + N_blanks))

    # Initialize vector to track call times
    call_times <- c()

    # Set initial parameters
    flip_rate_10_current <- runif(n = N_bits, min = 0.01, max = 0.025)
    flip_rate_01_current <- runif(n = N_bits, min = 0.01, max = 0.025)
    lum_noise_correlations_current <- runif(n = N_corrs, min = -0.05, max = 0.05)

    # Run simulation with initial parameters
    cat("\n\nRunning initial simulation with random bit-flip rates and luminance-noise correlations")
    sim_results <- run_simulation(
      codebook = codebook,
      summary_stats_genes = summary_stats_genes,
      summary_stats_blanks = summary_stats_blanks,
      flip_rate_10 = flip_rate_10_current,
      flip_rate_01 = flip_rate_01_current,
      max_correctable_Hamming_distance = 4,
      bit_lum_noise_correlations = lum_noise_correlations_current,
      blank_weight = max_blank_weight,
      maintain_gene_identity = FALSE,
      return_flip_rates = FALSE,
      verbose = FALSE
    )

    # Extract sim_mse of the initial simulation
    sim_mse_current <- sim_results$sim_summary["sim_mse"]
    cat("\nInitial mse:", round(sim_mse_current, 3))

    cat("\n\nRunning Markov chain steps: \n\nRunning initial steps ...")
    maintain_gene_identity <- FALSE
    return_flip_rates <- FALSE
    last_message <- 0
    last_step <- 0
    last_ctr <- 0
    FR_sample_range <- c(as.integer(FR_search*0.8):FR_search)
    progress_tracker <- as.integer(seq(1, n_steps, length.out = n_tracker_updates + 1))
    progress_tracker[n_tracker_updates] <- n_steps
    flip_rate_10_mean <- 1
    flip_rate_01_mean <- 1
    while (step < n_steps) {

      # Start timer
      start_time <- Sys.time()

      # Advance general counter
      ctr <- ctr + 1

      # Generate random step (... this is the Markov chain)
      flip_rate_10_next <- rnorm(n = N_bits, mean = flip_rate_10_current, sd = step_size_schedule[step + 1])
      flip_rate_01_next <- rnorm(n = N_bits, mean = flip_rate_01_current, sd = step_size_schedule[step + 1])
      lum_noise_correlations_next <- rnorm(N_corrs, mean = lum_noise_correlations_current, sd = step_size_schedule[step + 1])

      # Stochastically prune down
      retained_bits_idx <- sort(sample(1:N_bits, n_retained_bits[step + 1], replace = FALSE), decreasing = FALSE)
      retained_corrs_idx <- sort(sample(1:N_corrs, n_retained_corrs[step + 1], replace = FALSE), decreasing = FALSE)
      if (length(retained_bits_idx) > 0) {
        flip_rate_10_next[retained_bits_idx] <- flip_rate_10_current[retained_bits_idx]
        flip_rate_01_next[retained_bits_idx] <- flip_rate_01_current[retained_bits_idx]
      }
      if (length(retained_corrs_idx) > 0) {
        lum_noise_correlations_next[retained_corrs_idx] <- lum_noise_correlations_current[retained_corrs_idx]
      }

      # If burnin, enforce mean flip rates from flip-rate search
      # ... Note: In the MCMCSA Flip Rates plot, the clamped value in the middle step block likely won't match the tail of the
      #      value in the first step block, because this plot shows stipulated flip rate, while the middle block is clamped to the
      #      observed mean flip rate in the tail of the first block, and the stipulated and observed mean flip rates can be different.
      if (step > FR_search && step <= FR_search + burnin) {
        if (!is.na(flip_rate_10_mean) && !is.na(flip_rate_01_mean)) {
          flip_rate_10_next <- flip_rate_10_next - mean(flip_rate_10_next) + flip_rate_10_mean
          flip_rate_01_next <- flip_rate_01_next - mean(flip_rate_01_next) + flip_rate_01_mean
        }
      }

      # Enforce flip-rate bounds of 0 and 1 ("oob" = "out of bounds")
      oob_mask_10_next <- flip_rate_10_next < 0 | flip_rate_10_next > 1
      oob_mask_01_next <- flip_rate_01_next < 0 | flip_rate_01_next > 1
      if (step > FR_search && step <= FR_search + burnin) {
        if (!is.na(flip_rate_10_mean) && !is.na(flip_rate_01_mean)) {
          flip_rate_10_next[oob_mask_10_next] <- flip_rate_10_mean
          flip_rate_01_next[oob_mask_01_next] <- flip_rate_01_mean
        }
      } else {
        flip_rate_10_next[oob_mask_10_next] <- flip_rate_10_current[oob_mask_10_next]
        flip_rate_01_next[oob_mask_01_next] <- flip_rate_01_current[oob_mask_01_next]
      }
      # Enforce luminance noise correlation bounds of -1 and 1
      oob_mask_next <- lum_noise_correlations_next < -1 | lum_noise_correlations_next > 1
      lum_noise_correlations_next[oob_mask_next] <- lum_noise_correlations_current[oob_mask_next]

      # Run simulation with next parameters
      # ... this is step 1 of the Monte Carlo method: performing deterministic calculation of FPQC metrics and simulation-to-data fit
      sim_results <- run_simulation(
        codebook = codebook,
        summary_stats_genes = summary_stats_genes,
        summary_stats_blanks = summary_stats_blanks,
        flip_rate_10 = flip_rate_10_next,
        flip_rate_01 = flip_rate_01_next,
        max_correctable_Hamming_distance = 4,
        bit_lum_noise_correlations = lum_noise_correlations_next,
        blank_weight = blank_weight_schedule[step + 1],
        maintain_gene_identity = maintain_gene_identity,
        return_flip_rates = return_flip_rates,
        verbose = FALSE
      )

      # Extract log-likelihood of the next simulation
      sim_mse_next <- sim_results$sim_summary["sim_mse"]

      # Use simulated annealing to decide whether to accept or reject the proposed step

      # Calculate acceptance probability
      # ... idea: When mse decreases, probability of acceptance is 1; this formula
      #      controls how quickly the probability of acceptance decreases as the mse increases
      acceptance_prob <- min(1,exp(-(sim_mse_next - sim_mse_current)/temperature_schedule[step + 1]))

      # Accept or reject the proposed step
      ran_draw <- runif(n = 1, min = 0, max = 1)
      if (ran_draw < acceptance_prob || step == FR_search) {
        # Accept the new parameters
        # ... this is updating for the Markov chain
        flip_rate_10_current <- flip_rate_10_next
        flip_rate_01_current <- flip_rate_01_next
        lum_noise_correlations_current <- lum_noise_correlations_next
        # Update the sim_mse
        sim_mse_current <- sim_mse_next
        # Advance step
        step <- step + 1
        # Save results
        # ... this is step 2 of the Monte Carlo method: aggregate results
        sim_summaries[step,] <- sim_results$sim_summary
        PPV_genes[step,] <- sim_results$gene_summary$PPV
        if (step == 1) colnames(PPV_genes) <- rownames(sim_results$gene_summary)
        Counts_sorted_sim[step,] <- sim_results$counts_sorted$sim
        Counts_sorted_obs[step,] <- sim_results$counts_sorted$obs
        if (step == FR_sample_range[1] - 1 || step == FR_sample_range[length(FR_sample_range)] + 1) return_flip_rates <- !return_flip_rates
        # Grab mean flip rates
        if (step == FR_search) {
          # ... flip flag
          maintain_gene_identity <- TRUE
          # ... grab mean flip rates
          flip_rate_10_mean <- mean(sim_summaries[FR_sample_range, names(sim_results$sim_summary) == "flip_rate_10"], na.rm = TRUE)
          flip_rate_01_mean <- mean(sim_summaries[FR_sample_range, names(sim_results$sim_summary) == "flip_rate_01"], na.rm = TRUE)
          # ... reset parameters
          flip_rate_10_current <- runif(n = N_bits, min = 0.01, max = 0.025)
          flip_rate_01_current <- runif(n = N_bits, min = 0.01, max = 0.025)
          lum_noise_correlations_current <- runif(n = N_corrs, min = -0.05, max = 0.05)
        }
      }

      # End timer
      diff <- Sys.time() - start_time
      units(diff) <- "secs"
      call_times <- c(call_times, diff)

      # Update progress:
      if ((any(step == progress_tracker) || step == 10 || step == FR_search) && step != 1 && step != last_message) {
        # ... compute acceptance rate for this step
        acceptance_rate <- (step - last_step) / (ctr - last_ctr)
        last_step <- step
        last_ctr <- ctr
        if (acceptance_rate < 0.2) {
          temperature_schedule <- temperature_schedule * 1.1
          step_size_schedule <- step_size_schedule * 0.9
        } else if (acceptance_rate > 0.3) {
          temperature_schedule <- temperature_schedule * 0.9
          step_size_schedule <- step_size_schedule * 1.1
        }
        # ... report
        cat("\n\nStep", step, "/", n_steps, "complete")
        cat("\n\tCurrent mse:", round(sim_mse_current, 3))
        cat("\n\tAcceptance rate (this batch):", round(acceptance_rate, 3))
        cat("\n\tAcceptance rate (overall):", round(step / ctr, 3))
        cat("\n\tCalls:", ctr)
        cat("\n\tMean time (sec) per call:", round(mean(call_times), 3))
        if (step == FR_search) {
          cat("\n\nFlip-rate search complete, mean flip rates:",
              "\n\t1->0:", round(flip_rate_10_mean, 3),
              "\n\t0->1:", round(flip_rate_01_mean, 3))
        }
        last_message <- step
      }

    }

    colnames(sim_summaries) <- names(sim_results$sim_summary)

    # Make PPV plot to evaluate run
    resamples_rows <- c((FR_search + burnin):(FR_search + burnin + resamples))
    PPV_plot_info <- plot_PPV_cutoff(PPV_genes[resamples_rows, ])
    PPV_summary <- PPV_plot_info$summary_df
    PPV_plot <- PPV_plot_info$PPV_plot
    print(PPV_plot)

    # Make sorted count plot to evaluate simulation
    counts_sorted_sim2 <- Counts_sorted_sim[resamples_rows,]
    counts_sorted_sim_ <- apply(counts_sorted_sim2, 2, median)
    counts_sorted_obs_ <- Counts_sorted_obs[1,]
    counts_sorted_obs_ <- unlist(as.vector(counts_sorted_obs_))
    sorted_count_plot <- plot_counts_sorted(counts_sorted_sim_, counts_sorted_obs_, N_genes)
    print(sorted_count_plot)

    # Plot steps to evaluate the resampling
    plot_MCMC_steps(sim_summaries, burnin, FR_search)

    # Bundle results
    test_results <- list(
      sim_summaries = sim_summaries,
      PPV = PPV_genes,
      PPV_summary = PPV_summary,
      Counts_sorted_sim = Counts_sorted_sim,
      Counts_sorted_obs = Counts_sorted_obs,
      PPV_plot = PPV_plot,
      sorted_count_plot = sorted_count_plot
    )

    # Save to .rds file
    # MCMCSA_test <- readRDS("FPQC_run_m1.rds")
    saveRDS(test_results, file = paste0("FPQC_run_", mouse_id, ".rds"))
    write.csv(test_results$PPV_summary, file = paste0("FPQC_run_", mouse_id, "_PPV_summary.csv"), row.names = FALSE)

    return(test_results)

  }

# Plot run QC (PPV cutoff)
plot_PPV_cutoff <- function(
    PPV_resamples,
    CI_range = 0.95,
    min_PPV = 0.8
  ) {
    
    CI_range_low <- (1 - CI_range)/2
    CI_range_high <- CI_range + CI_range_low
    
    # Compute column medians
    #medians <- apply(PPV_resamples, 2, median)
    # Order column names by decreasing median
    #ordered_cols <- names(sort(medians, decreasing = TRUE))
    lower_q <- apply(PPV_resamples, 2, function(x) quantile(x, CI_range_low, na.rm = TRUE))
    # Order columns by decreasing lower quantile
    ordered_cols <- names(sort(lower_q, decreasing = TRUE))
    
    # Reorder columns in PPV_resamples
    PPV_resamples_ordered <- PPV_resamples[ , ordered_cols]
    col_range <- sample(1:ncol(PPV_resamples_ordered), size = ncol(PPV_resamples_ordered), replace = FALSE)
    col_range <- sort(col_range, decreasing = FALSE)
    PPV_resamples_long <- stack(PPV_resamples_ordered)
    colnames(PPV_resamples_long) <- c("value", "variable")
    
    # Summarize median and 95% confidence interval
    
    summary_df <- PPV_resamples_long %>%
      group_by(variable) %>%
      summarize(
        median = median(value, na.rm = TRUE),
        lower = quantile(value, CI_range_low, na.rm = TRUE),
        upper = quantile(value, CI_range_high, na.rm = TRUE)
      )
    summary_df$above_cutoff <- rep("bad", nrow(summary_df))
    summary_df$above_cutoff[summary_df$lower > min_PPV] <- "good"
    
    # Plot
    PPV_plot <- ggplot(summary_df, aes(x = variable, y = median)) +
      geom_point(aes(color = above_cutoff), size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper, color = above_cutoff), width = 0.2) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 18),      # Axis titles
        axis.text = element_text(size = 12),       # Axis tick labels
        plot.title = element_text(size = 24, hjust = 0.5),  # Title
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
      scale_color_manual(values = c("good" = "blue", "bad" = "black")) +
      guides(color = "none") +
      labs(
        y = paste0("Estimated PPV: Median with ", CI_range, " CI"),
        x = paste0("Gene barcodes by lower ", CI_range_low, " quantile of PPV confidence interval"),
        title = "Estimated Precision (PPV) by Gene Barcode")
    # ^ Save as 1600 x 925
    
    return(
      list(PPV_plot = PPV_plot, summary_df = summary_df)
    )
    
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

# 6. Run and analyze simulations #######################################################################################

# MCMCSA_test <- run_MCMCSA(
#   codebook = codebook,
#   mouse_id = "m1", # m1, m2, m3, etc...
#   resamples = 1000,
#   burnin = 1000,
#   FR_search = 1000
#   )

remake_plots <- function() {
  files <- list.files(pattern = "\\.rds$", full.names = TRUE)
  rds_list <- lapply(files, readRDS)
  burnin <- 1000
  resamples <- 1000
  FR_search <- 1000
  resamples_rows <- c((FR_search + burnin):(FR_search + burnin + resamples))
  for (i in seq_along(rds_list)) {
    sim <- rds_list[[i]]
    PPV_plot_info <- plot_PPV_cutoff(sim[["PPV"]][resamples_rows,])
    PPV_summary <- PPV_plot_info$summary_df
    PPV_plot <- PPV_plot_info$PPV_plot
    print(PPV_plot)
    N_genes <- nrow(sim[["PPV_summary"]])
    counts_sorted_sim_ <- apply(sim[["Counts_sorted_sim"]][resamples_rows,], 2, median)
    counts_sorted_obs_ <- sim[["Counts_sorted_obs"]][1,]
    counts_sorted_obs_ <- unlist(as.vector(counts_sorted_obs_))
    sorted_count_plot <- plot_counts_sorted(counts_sorted_sim_, counts_sorted_obs_, N_genes)
    print(sorted_count_plot)
  }
}
#remake_plots()


# So far, this is the best settings (see defaults); next, try starting with lower flip-rate range (0>1 does better and starts at .15, 1>0 does worse and starts at .3)
# MSE log isn't having a long stable run at the end; adjust the cooling and step-size schedule to a flatter curve

#write.csv(MCMCSA_test$PPV, file = "MCMCSA_test_PPV.csv", row.names = FALSE)
#write.csv(MCMCSA_test$sim_summaries, file = "MCMCSA_test_sim_summaries.csv", row.names = FALSE)

#write.csv(sim_results, file = "sim_results.csv", row.names = FALSE)
#sim_results <- read.csv("sim_results.csv")






