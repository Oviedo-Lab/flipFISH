
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

