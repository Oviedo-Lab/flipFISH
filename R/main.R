
# 1. Setup #############################################################################################################

# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2026

#' @useDynLib flipFISH, .registration = TRUE
#' @import Rcpp
#' @import RcppEigen 
#' @import ggplot2
NULL

.onLoad <- function(libname, pkgname) {}

#' Run misread QC analysis on summary stats data from spatial transcriptomics experiment
#'
#' This function takes summary statistics from a FISH-based spatial transcriptomics experiment and the barcode codebook (including blanks labelled with "Blank") and runs L-BFGS (via nlopt) to estimate the bit-flip rates and correlations, as well as the expected read, error-corrected, and true counts for each barcode. The aim is to compute both the "confidence ratio" (CR) and positive predictive value (PPV) for each barcode.
#'
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts", must have barcode names as row names
#' @param codebook Codebook with row names as barcodes and columns as bits, must have barcode names as row names
#' @param max_fr Maximum flip rate allowed during optimization, default 1.0
#' @param max_corr Absolute bound on bit-flip correlation parameters, default 1.0
#' @param initial_corr_scale Assumed mean correlation between bit-flips, used to make crude estimate of true barcode counts from observed counts, default 0.1
#' @param ah_hoc_rescale Additional rescaling of observed counts when estimating true counts, default is 1.1
#' @param n_forks Number of parallel forks to use for expected-count computation, default is 1 (must be 1 on Windows)
#' @param max_flips When computing expected counts per barcode, ignore misreads larger than this hamming distance, default is 0, interpreted as no limit
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every report_freq accepted calls, default 1
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 1000
#' @param ftol_rel Relative function value tolerance for L-BFGS convergence, default 1e-8
#' @param xtol_rel Relative parameter tolerance for L-BFGS convergence, default 1e-6
#' @param max_correctable_Hamming_distance Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries
#' @return List giving \code{STdata}, a dataframe of the ST summary data, \code{fliprates}, a 1-row matrix of the optimal flip rates and bit-flip correlations, \code{erc}, \code{ecc}, and \code{etc}, 1-row matrices of the expected read, error-corrected, and true counts at the optimum, \code{CR} and \code{PPV}, 1-row matrices of the confidence ratio and PPV at the optimum, \code{msle}, the final objective value, and \code{fliprates_summary} and \code{bc_summary}, summary tables (point estimates; lower == upper == mean since L-BFGS returns a single solution).
#' @export
misreadQC <- function(
    STdata,
    codebook,
    max_fr                           = 1.0,
    max_corr                         = 1.0,
    initial_corr_scale               = 0.1, 
    ah_hoc_rescale                   = 1.1, 
    n_resamples                      = 100,
    n_forks                          = 1,
    max_flips                        = 0,
    report_freq                      = 1,
    maxeval                          = 1000,
    ftol_rel                         = 1e-8,
    xtol_rel                         = 1e-6,
    max_correctable_Hamming_distance = NULL
  ) {
    cat("\nRunning misread QC with L-BFGS (nlopt)")
    cat("\nMax flip rate:", max_fr)
    cat("\nMax evaluations:", maxeval)
    
    # Confirm forking is possible and check number of cores
    if (!(Sys.info()["sysname"] == "Darwin" || Sys.info()["sysname"] == "Linux")) {
      if (n_forks > 1) {
        cat("\nForking not available on Windows, setting n_forks to 1")
        n_forks <- 1
      }
    } else if (n_forks > parallel::detectCores()) {
      cat("\nn_forks exceeds available cores, setting n_forks to", parallel::detectCores())
      n_forks <- parallel::detectCores()
    } else {
      cat("\nNumber of forks to use:", n_forks)
    }
    
    # Prep codebook and STdata
    # ... get species names
    species_names <- rownames(STdata)
    if (is.null(species_names)) stop("STdata must have row names as species names")
    # ... align rows of codebook to STdata
    if (is.null(rownames(codebook))) stop("codebook must have row names as species names")
    if (!all(species_names %in% rownames(codebook))) stop("All species in STdata must be present in codebook")
    codebook <- codebook[species_names,]
    # ... make blank mask
    blank_mask <- grepl("Blank", species_names, ignore.case = FALSE)
    if (sum(blank_mask) == 0) stop("No blanks found in STdata row names, make sure blank species have 'Blank' in their names")
    if (sum(blank_mask) == length(species_names)) stop("All species are blanks. Make sure non-blank species do not have 'Blank' in their names")
    # ... sort species by decreasing rates, with genes first and blanks second 
    gene_rate_order <- order(STdata$rates[!blank_mask], decreasing = TRUE)
    blank_rate_order <- order(STdata$rates[blank_mask], decreasing = TRUE)
    # ... remake STdata and codebook with this order
    STdata <- rbind(
      STdata[!blank_mask,][gene_rate_order,],
      STdata[blank_mask,][blank_rate_order,]
    )
    codebook <- rbind(
      codebook[!blank_mask,][gene_rate_order,],
      codebook[blank_mask,][blank_rate_order,]
    )
    
    # Check bit size 
    if (ncol(codebook) > 64) stop("Codebook has more than the max-allowed 64 bits.")
    
    # Check and set max_correctable_Hamming_distance
    codebook_distances <- unique_Hamming_cb(as.matrix(codebook))
    if (is.null(max_correctable_Hamming_distance)) {
      max_correctable_Hamming_distance <- min(codebook_distances) - 1
    } else if (max_correctable_Hamming_distance >= min(codebook_distances)) {
      stop(paste0("max_correctable_Hamming_distance must be less than the minimum Hamming distance between codebook entries (", min(codebook_distances), ")"))
    }
    
    # Run misread QC algorithm with L-BFGS
    qc <- mQC(
      as.matrix(STdata),
      as.matrix(codebook),
      max_correctable_Hamming_distance,
      max_fr,
      max_corr,
      initial_corr_scale,
      ah_hoc_rescale,
      as.integer(n_resamples), 
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      maxeval,
      ftol_rel,
      xtol_rel
    )
    
    # Make summary stats from qc results
    cat("\nRunning summary stats on QC results")
    sum_names <- c("mean", "lower", "upper")
    qc_names <- names(qc)
    bc_names <- qc_names[qc_names != "STdata" & qc_names != "fliprates"]
    bc_sum_names <- c()
    for (n in bc_names) {bc_sum_names <- c(bc_sum_names, paste0(n, "_", sum_names))}
    fr <- matrix(NA, nrow = ncol(qc$fliprates), ncol = length(sum_names))
    bc <- matrix(NA, nrow = ncol(qc$ecc), ncol = length(sum_names) * length(bc_names))
    colnames(fr) <- sum_names 
    colnames(bc) <- bc_sum_names
    rownames(bc) <- qc$STdata$species
    N_bits <- ncol(codebook)
    fr_names <- paste0("rate10_bit", seq_len(N_bits))
    fr_names <- c(fr_names, paste0("rate01_bit", seq_len(N_bits)))
    fr_names <- c(fr_names, paste0("corr_", seq_len(ncol(qc$fliprates) - 2*N_bits)))
    rownames(fr) <- fr_names
    for (p in qc_names) {
      if (p == "STdata" || p == "msle") next
      if (n_resamples == 1) {
        p_means <- qc[[p]]
        ci <- rbind(p_means, p_means)
      } else {
        p_means <- colMeans(qc[[p]])
        ci <- apply(qc[[p]], 2, quantile, probs = c(0.025, 0.975))
      }
      if (p == "fliprates") {
        fr[,"mean"] <- p_means
        fr[,"lower"] <- ci[1,]
        fr[,"upper"] <- ci[2,]
      } else {
        bc[,paste0(p, "_mean")] <- p_means
        bc[,paste0(p, "_lower")] <- ci[1,]
        bc[,paste0(p, "_upper")] <- ci[2,]
      }
    }
    qc[["fliprates_summary"]] <- fr
    qc[["bc_summary"]]        <- bc
    
    rate10_mean <- mean(fr[grepl("rate10", rownames(fr)), "mean"])
    rate01_mean <- mean(fr[grepl("rate01", rownames(fr)), "mean"])
    rate10_lower <- mean(fr[grepl("rate10", rownames(fr)), "lower"])
    rate10_upper <- mean(fr[grepl("rate10", rownames(fr)), "upper"])
    rate01_lower <- mean(fr[grepl("rate01", rownames(fr)), "lower"])
    rate01_upper <- mean(fr[grepl("rate01", rownames(fr)), "upper"])
    cat("\nEstimated flip rates (mean, 95% CI):")
    cat("\n1>0: ", round(rate10_mean, 4), " (", round(rate10_lower, 4), "-", round(rate10_upper, 4), ")", sep = "")
    cat("\n0>1: ", round(rate01_mean, 4), " (", round(rate01_lower, 4), "-", round(rate01_upper, 4), ")\n", sep = "")
    
    return(qc)
    
  }

#' Plot estimated PPV by barcode from misread QC results
#'
#' This function takes the results from the misreadQC function and makes a plot of the estimated positive predictive value (PPV) for each barcode, with error bars showing the 95% confidence intervals. Barcodes are sorted by decreasing lower bound of the PPV confidence interval, and a dashed red line indicates the minimum PPV cutoff for "good" barcodes. Barcodes whose entire confidence interval is above the cutoff are colored blue, while those that are not are colored black. 
#'
#' @name plot.PPV
#' @rdname plot-PPV
#' @usage plot.PPV(
#'  qc,
#'  min_PPV = 0.8
#' )
#' @param qc List, results from misreadQC function
#' @param min_PPV Numeric, minimum 95% CI lower-bound PPV cutoff for "good" barcodes, defaults to 0.8
#' @return ggplot object showing estimated PPV for each barcode with confidence intervals and cutoff
#' @export
plot.PPV <- function(
    qc, # List, results from misreadQC function
    min_PPV = 0.8
  ) {
    
    PPV_mean <- qc$bc_summary[,"PPV_mean"]
    PPV_ci_lower <- qc$bc_summary[,"PPV_lower"]
    PPV_ci_upper <- qc$bc_summary[,"PPV_upper"]
    
    # Sort by decreasing values
    sorted_idx <- order(PPV_ci_lower, decreasing = TRUE)
    PPV_mean <- PPV_mean[sorted_idx]
    PPV_ci_lower <- PPV_ci_lower[sorted_idx]
    PPV_ci_upper <- PPV_ci_upper[sorted_idx]
    PPV_cutoff <- rep("bad", length(PPV_mean))
    PPV_cutoff[PPV_ci_lower > min_PPV] <- "good"
    df <- data.frame(
      x = seq_along(PPV_mean), 
      PPV = PPV_mean, 
      lower = PPV_ci_lower,
      upper = PPV_ci_upper,
      above_cutoff = PPV_cutoff
      )
    
    # Plot
    PPV_plot <- ggplot(df, aes(x = x, y = PPV)) +
      geom_point(aes(color = above_cutoff), size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper, color = above_cutoff), width = 0.2) +
      geom_hline(yintercept = min_PPV, linetype = "dashed", color = "red") +
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
    
    return(PPV_plot)
    
  }

#' Plot predicted vs observed counts from misread QC results
#' 
#' This function takes the results from the misreadQC function and makes a plot comparing the predicted error-corrected counts (ecc) to the observed counts for each barcode. Both predicted and observed counts are plotted on a log scale for better visibility, with error bars showing the 95% confidence intervals for the predicted counts. Barcodes are sorted by decreasing observed count, with gene barcodes shown first and blank barcodes shown second. A dashed vertical line indicates the separation between gene and blank barcodes.
#' 
#' @name plot.counts 
#' @rdname plot-counts
#' @usage plot.counts(qc)
#' @param qc List, results from misreadQC function
#' @return ggplot object showing predicted vs observed counts for each barcode with confidence intervals
#' @export
plot.counts <- function(
    qc
  ) {
    
    bc_type_colors <- c(
      "Gene (pred)" = "skyblue1", "Blank (pred)" = "gray",
      "Gene (obs)" = "skyblue4", "Blank (obs)" = "gray20"
    )
    
    count_obs <- qc[["STdata"]]$count_observed
    count_pred <- qc[["bc_summary"]][,"ecc_mean"] 
    count_lower <- qc[["bc_summary"]][,"ecc_lower"]
    count_upper <- qc[["bc_summary"]][,"ecc_upper"]
    
    blank_mask <- grepl("Blank", qc[["STdata"]]$species)
    bc_type <- rep("Gene (pred)", length(count_pred))
    bc_type[blank_mask] <- "Blank (pred)"
    bc_type_obs <- rep("Gene (obs)", length(count_obs))
    bc_type_obs[blank_mask] <- "Blank (obs)"
    
    df <- data.frame(
      index = c(seq_along(count_pred), seq_along(count_obs)), 
      Count = c(count_pred, count_obs),
      bc_type = c(bc_type, bc_type_obs),
      size = c(rep(1.5, length(count_pred)), rep(1.5, length(count_obs)))
    )
    df$Count[df$Count == 0] <- 1
    df$bc_type <- factor(
      df$bc_type, 
      levels = c(
        "Gene (pred)", "Gene (obs)", 
        "Blank (pred)", "Blank (obs)"
      )
    )
    
    df_pred <- data.frame(
      index = seq_along(count_pred),
      Count = count_pred,
      lower = count_lower,
      upper = count_upper,
      bc_type = bc_type
    )
    df_pred$Count[df_pred$Count == 0] <- 1
    df_pred$bc_type <- factor(
      df_pred$bc_type, 
      levels = c(
        "Gene (pred)", "Gene (obs)", 
        "Blank (pred)", "Blank (obs)"
      )
    )
    
    N_genes <- sum(!blank_mask)
    counts_sorted_plot <- ggplot(df) +
      geom_errorbar(data = df_pred, aes(x = index, ymin = lower, ymax = upper, color = bc_type), width = 0.2, alpha = 0.5) +
      geom_point(size = df$size, aes(x = index, y = Count, color = bc_type)) +
      geom_vline(
        xintercept = mean(c(N_genes, N_genes+1)), 
        color = "black", 
        linetype = "dashed", 
        linewidth = 0.5
      ) +
      scale_y_log10() +  # Log scale for better visibility
      scale_color_manual(values = bc_type_colors) +
      labs(title = "Predicted vs Observed Counts", x = "Barcodes sorted by observed count", y = "Spot count", color = "Count type") +
      theme_minimal() +
      theme(
        legend.position = "right"
      )
    
    return(counts_sorted_plot)
    
  }

#' Plot estimated flip rate distributions from misread QC results
#' 
#' This function takes the results from the misreadQC function and makes a plot showing the distributions of the estimated bit-flip rates for each bit, separated by flip type (1>0 vs 0>1). The plot uses violin plots to show the distribution of flip rates across the second half of iterations of the MCMCSA walk, with separate facets for each flip type.
#' 
#' @name plot.fr
#' @rdname plot-fr
#' @usage plot.fr(qc)
#' @param qc List, results from misreadQC function
#' @return ggplot object showing distributions of estimated flip rates for each bit and flip type
#' @export
plot.fr <- function(
    qc
  ) {
    fr_names <- rownames(qc$fliprates_summary)
    mask10 <- grepl("rate10", fr_names)
    mask01 <- grepl("rate01", fr_names)
    n_samples <- nrow(qc$fliprates)
    N_bits <- sum(mask10)
    if (sum(mask01) != N_bits) stop("Number of rate10 and rate01 entries in fliprates_summary must be the same")
    fr <- qc$fliprates
    fr <- matrix(NA, nrow = 2*n_samples*N_bits, ncol = 3)
    colnames(fr) <- c("value", "bit", "type")
    fr <- as.data.frame(fr)
    for (i in seq_len(N_bits)) {
      idx10 <- (i-1)*n_samples + seq_len(n_samples)
      idx01 <- (N_bits + i-1)*n_samples + seq_len(n_samples)
      fr[idx10, "value"] <- qc$fliprates[,paste0("rate10_bit", i) == fr_names]
      fr[idx10, "bit"] <- i
      fr[idx10, "type"] <- "1>0"
      fr[idx01, "value"] <- qc$fliprates[,paste0("rate01_bit", i) == fr_names]
      fr[idx01, "bit"] <- i
      fr[idx01, "type"] <- "0>1"
    }
    fr$bit <- as.factor(fr$bit)
    plt <- ggplot(fr, aes(bit, value, fill = type)) +
      geom_violin() +
      labs(title = "Estimated Flip Rate Distributions", x = "Bit", y = "Flip Rate", color = "Flip Type") +
      theme_minimal() +
      facet_grid(type ~ .)
    return(plt)
  }
