
# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2026

#' @useDynLib flipFISH, .registration = TRUE
#' @import Rcpp
#' @import RcppEigen 
#' @import ggplot2
NULL

.onLoad <- function(libname, pkgname) {}

# Helper function to check STdata and codebook 
check_data <- function(
    STdata, 
    codebook,
    max_correctable_Hamming_distance
  ) {
    # Check bit size 
    if (ncol(codebook) > 64)                         stop("Codebook has more than the max-allowed 64 bits.")
   
    # Prep codebook and STdata
    # ... get species names
    species_names <- rownames(STdata)
    if (is.null(species_names))                      stop("STdata must have row names as barcode species names")
    # ... align rows of codebook to STdata
    if (is.null(rownames(codebook)))                 stop("codebook must have row names as species names")
    if (!all(species_names %in% rownames(codebook))) stop("All species in STdata must be present in codebook")
    codebook      <- codebook[species_names,]
    # ... make blank mask
    blank_mask    <- grepl("Blank", species_names, ignore.case = FALSE)
    if (sum(blank_mask) == 0)                        stop("No blanks found in STdata row names, make sure blank species have 'Blank' in their names")
    if (sum(blank_mask) == length(species_names))    stop("All species are blanks. Make sure non-blank species do not have 'Blank' in their names")
    # ... sort species by decreasing rates, with genes first and blanks second 
    gene_rate_order  <- order(STdata$rates[!blank_mask], decreasing = TRUE)
    blank_rate_order <- order(STdata$rates[blank_mask],  decreasing = TRUE)
    # ... remake STdata and codebook with this order
    STdata   <- rbind(
      STdata[!blank_mask,][gene_rate_order, ],
      STdata[ blank_mask,][blank_rate_order,]
    )
    codebook <- rbind(
      codebook[!blank_mask,][gene_rate_order, ],
      codebook[ blank_mask,][blank_rate_order,]
    )
    
    # Check and set max_correctable_Hamming_distance
    codebook_distances <- unique_Hamming_cb(as.matrix(codebook))
    if (is.null(max_correctable_Hamming_distance)) {
      max_correctable_Hamming_distance <- min(codebook_distances) - 1
    } else if (max_correctable_Hamming_distance >= min(codebook_distances)) {
      stop(paste0("max_correctable_Hamming_distance must be less than the minimum Hamming distance between codebook entries (", min(codebook_distances), ")"))
    }
    
    return(
      list(
        STdata                           = STdata,
        codebook                         = codebook,
        max_correctable_Hamming_distance = max_correctable_Hamming_distance
      )
    )
  }

# Helper function to check forks 
check_forks <- function(
    n_forks
  ) {
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
    return(n_forks)
  }

#' Run misread QC analysis on summary stats data from spatial transcriptomics experiment
#'
#' This function takes summary statistics from a FISH-based spatial transcriptomics experiment and the barcode codebook (including blanks labelled with "Blank") and runs L-BFGS (via nlopt) to make a best-fit estimate of the bit-flip rates and correlations. The estimation (i.e., L-BFGS optimization) uses an analytic conditional probability model to compute, for each barcode \emph{b}, values expected based on bit-flip rates and bit-flip correlations, for: \itemize{
#'   \item \strong{Read Count}: The number of spots (the "count") labelled \emph{b} before error correction.
#'   \item \strong{Corrected Count}: The number of spots (the "count") labelled \emph{b} after error correction.
#'   \item \strong{Hit Count}: The number of spots (the "count") labelled \emph{b} after error correction which are in fact \emph{b}, i.e., the number of reads error-corrected to \emph{b} which are correct, aka a "hit". 
#' }
#'  The "best fit" is defined as least mean squared log error, mean squared log error being computed by comparing the analytically implied expected corrected counts to the observed corrected counts included in the summary statistics. Additionally, both the expected "confidence ratio" (CR) and expected positive predictive value (PPV) are computed, from these implied values, for each barcode. The CR is a quality control metric proposed by the original developers of MERFISH (e.g., see DOI 10.1016/bs.mie.2016.03.020) defined per barcode as the read count over the corrected count, while PPV is a oft-used metric defined as the hit count over the corrected count. That is, PPV tells us, for each barcode \emph{b}, the expected percentage of spots decoded as \emph{b} which are in fact \emph{b}. 
#'
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}. 
#' @param n_forks Number of process forks to use for expected-count computation (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param max_flips When analytically computing expected corrected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between six and ten is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing. 
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
#' @param max_correctable_Hamming_distance Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries.
#' @return A list giving:\itemize{
#'    \item \code{STdata}: a dataframe giving the summary data from the ST run used in the estimation.
#'    \item \code{fliprates}: a labeled vector giving the estimated (i.e., best fit) flip rates and bit-flip correlations from the optimization. 
#'    \item \code{erctc_plus}: a data frame giving, for each barcode, the read, corrected, and hit counts as well as \code{CR} and \code{PPV} values implied by the analytic conditional probability model, given the values in \code{fliprates}. Column names are: "erc", "ecc", "etc", "CR", and "PPV".
#'    \item \code{msle}: A numeric value giving the minimal mean squared log error found by the L-BFGS optimization of the flip rates and bit-flip correlations.
#' } 
#' @export
misread.qc <- function(
    STdata,
    codebook,
    n_forks                          = 1,
    max_flips                        = 0,
    report_freq                      = 10,
    maxeval                          = 500,
    max_correctable_Hamming_distance = NULL
  ) {
    cat("\nRunning misread QC with L-BFGS (nlopt)")
    cat("\nMax evaluations:", maxeval)
    
    # Confirm forking is possible and check number of cores
    n_forks    <- check_forks(n_forks)
    
    # Check bit size, prep codebook and STdata, and set max_correctable_Hamming_distance
    data_check <- check_data(STdata, codebook, max_correctable_Hamming_distance)
    
    # Run misread QC algorithm with L-BFGS
    qc <- mQC(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(data_check$max_correctable_Hamming_distance),
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      as.integer(maxeval),
      list()
    )
    
    # Annotate result output
    cat("\nBuilding summary tables")
    qc$erctc_plus            <- as.data.frame(qc$erctc_plus)
    row.names(qc$erctc_plus) <- qc$STdata$species
    N_bits                   <- ncol(codebook)
    names(qc$fliprates)      <- c(
      paste0("rate10_bit", seq_len(N_bits)),
      paste0("rate01_bit", seq_len(N_bits)),
      paste0("corr_",      seq_len(length(qc$fliprates) - 2*N_bits))
    )
    
    # Report flip-rate means
    rate10_mean <- mean(qc$fliprates[grepl("rate10", names(qc$fliprates))])
    rate01_mean <- mean(qc$fliprates[grepl("rate01", names(qc$fliprates))])
    cat("\nEstimated flip rates:")
    cat("\n1>0:", round(rate10_mean, 4))
    cat("\n0>1:", round(rate01_mean, 4))
    
    return(qc)
    
  }

#' Benchmark \code{misread.qc} function with dichotomized-Gaussian simulations
#' 
#' This function takes the same summary statistics (\code{STdata}) and barcode codebook (\code{codebook}) as \code{misread.qc} and runs dichotomized-Gaussian simulations with stipulated bit-flip rates and bit-flip correlations in order to estimate how well the L-BFGS algorithm recovers the bit-flip rates and bit-flip correlations for the given data set and codebook. 
#' 
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}. 
#' @param n_sims Number of dichotomized-Gaussian simulations to run. The default is 100. 
#' @param n_forks Number of process forks to use for expected-count computation (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param max_flips When analytically computing expected corrected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between six and ten is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing. 
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
#' @param max_correctable_Hamming_distance Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries.
#' @return A list giving:\itemize{
#'    \item \code{fr_est}: A matrix giving the flip rates and bit-flip correlations estimated by \code{misread.qc} (columns), for each set of observed counts generated by dichotomized-Gaussian simulation (rows).
#'    \item \code{fr_stipulated}: A vector giving the stipulated flip rates and bit-flip correlations used to generate the dichotomized-Gaussian simulations. All simulations use the same stipulated values.
#'    \item \code{PPV_est}: A matrix giving positive predictive value (PPV) expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows). 
#'    \item \code{PPV_expected}: A vector giving the PPV values expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}. 
#'    \item \code{ecc_est}: A matrix giving the corrected count expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows).  
#'    \item \code{ecc_expected}: A vector giving the corrected count expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}. 
#'    \item \code{sim_counts}: A matrix giving the actual simulated count per barcode (columns) for each simulation (rows). 
#'    }
dichot.guass.benchmark <- function(
    STdata,
    codebook,
    n_sims                           = 100,
    n_forks                          = 1,
    max_flips                        = 0,
    report_freq                      = 10,
    maxeval                          = 500,
    max_correctable_Hamming_distance = NULL
  ) {
    cat("\nBenchmarking misread QC with dichotomized-Gaussian simulation")
    
    # Confirm forking is possible and check number of cores
    n_forks    <- check_forks(n_forks)
    
    # Check bit size, prep codebook and STdata, and set max_correctable_Hamming_distance
    data_check <- check_data(STdata, codebook, max_correctable_Hamming_distance)
    
    # Run misread QC algorithm with L-BFGS
    resids <- test_fr_recovery(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(n_sims), 
      as.integer(data_check$max_correctable_Hamming_distance),
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      as.integer(maxeval),
      list()
    )
    
    return(resids)
    
  }

#' Plot estimated PPV by barcode from \code{misread.qc} results
#'
#' This function takes the results from the misread.qc function and makes a plot of the estimated positive predictive value (PPV) for each barcode. Barcodes are sorted by PPV value, in decreasing order, and a dashed red line indicates the minimum PPV cutoff for "good" barcodes. Barcodes above the cutoff are colored blue, while those that are not are colored black. 
#'
#' @name plot.PPV
#' @rdname plot-PPV
#' @usage plot.PPV(
#'  qc,
#'  min_PPV = 0.8
#' )
#' @param qc List, results from misread.qc function.
#' @param min_PPV Numeric, PPV cutoff for "good" barcodes, defaults to 0.8.
#' @return ggplot object showing estimated PPV for each barcode with cutoff.
#' @export
plot.PPV <- function(
    qc,
    min_PPV = 0.8
  ) {
    
    # Grab data
    PPV <- qc$erctc_plus$PPV
    
    # Sort by decreasing values
    above_cutoff <- rep("Bad", length(PPV))
    above_cutoff[PPV > min_PPV] <- "Good"
    df <- data.frame(
      x            = seq_along(PPV), 
      PPV          = PPV[order(PPV, decreasing = TRUE)], 
      above_cutoff = above_cutoff
      )
    
    # Plot
    PPV_plot <- ggplot(df, aes(x = x, y = PPV)) +
      geom_point(aes(color = above_cutoff), size = 3) +
      geom_hline(yintercept = min_PPV, linetype = "dashed", color = "red") +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank()) +
      scale_color_manual(values = c("Good" = "blue", "Bad" = "black")) +
      guides(color = "none") +
      labs(
        y     = paste0("Estimated PPV"),
        x     = paste0("Gene barcodes by PPV value"),
        title = "Estimated Precision (PPV) by Gene Barcode")
    
    return(PPV_plot)
    
  }

#' Plot predicted vs observed counts from \code{misread.qc} results
#' 
#' This function takes the results from the misread.qc function and makes a plot comparing the predicted error-corrected counts (ecc) to the observed counts for each barcode. Both predicted and observed counts are plotted on a log scale for better visibility. Barcodes are sorted by decreasing observed count, with gene barcodes shown first and blank barcodes shown second. A dashed vertical line indicates the separation between gene and blank barcodes.
#' 
#' @name plot.counts 
#' @rdname plot-counts
#' @usage plot.counts(qc)
#' @param qc List, results from misread.qc function.
#' @return ggplot object showing predicted vs observed counts for each barcode.
#' @export
plot.counts <- function(
    qc
  ) {
    
    # Set colors
    bc_type_colors <- c(
      "Gene (pred)" = "skyblue1", "Blank (pred)" = "gray",
      "Gene (obs)"  = "skyblue4", "Blank (obs)"  = "gray20"
    )
    
    # Grab data
    count_obs   <- qc$STdata$count_observed
    count_pred  <- qc$erctc_plus$ecc
    
    # Mask data
    blank_mask              <- grepl("Blank", qc[["STdata"]]$species)
    bc_type                 <- rep("Gene (pred)", length(count_pred))
    bc_type[blank_mask]     <- "Blank (pred)"
    bc_type_obs             <- rep("Gene (obs)", length(count_obs))
    bc_type_obs[blank_mask] <- "Blank (obs)"
    
    # Prepare data frame for plotting
    df <- data.frame(
      index   = c(seq_along(count_pred), seq_along(count_obs)), 
      Count   = c(count_pred, count_obs),
      bc_type = c(bc_type, bc_type_obs),
      size    = c(rep(1.5, length(count_pred)), rep(1.5, length(count_obs)))
    )
    df$Count[df$Count == 0] <- 1
    df$bc_type <- factor(
      df$bc_type, 
      levels = c(
        "Gene (pred)", "Gene (obs)", 
        "Blank (pred)", "Blank (obs)"
      )
    )
    
    # Prepare data frame for plotting
    df_pred <- data.frame(
      index   = seq_along(count_pred),
      Count   = count_pred,
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
    
    # Make plot
    N_genes <- sum(!blank_mask)
    counts_sorted_plot <- ggplot(df) +
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
      theme(legend.position = "right")
    
    return(counts_sorted_plot)
    
  }

#' Plot estimated flip rates from \code{misread.qc} results
#' 
#' This function takes the results from the misread.qc function and makes a plot showing the estimated bit-flip rates for each bit, separated by flip type (1>0 vs 0>1). 
#' 
#' @name plot.fr
#' @rdname plot-fr
#' @usage plot.fr(qc)
#' @param qc List, results from misread.qc function.
#' @return ggplot object showing distributions of estimated flip rates for each bit and flip type.
#' @export
plot.fr <- function(
    qc
  ) {
    # Make masks
    fr_names  <- rownames(qc$fliprates_summary)
    mask10    <- grepl("rate10", fr_names)
    mask01    <- grepl("rate01", fr_names)
    n_samples <- nrow(qc$fliprates)
    N_bits    <- sum(mask10)
    if (sum(mask01) != N_bits) stop("Number of rate10 and rate01 entries in fliprates_summary must be the same")
    # Grab data
    fr           <- matrix(NA, nrow = 2*n_samples*N_bits, ncol = 3)
    colnames(fr) <- c("value", "bit", "type")
    fr           <- as.data.frame(fr)
    for (i in seq_len(N_bits)) {
      # Indexing stipulated by the misread.qc function
      idx10              <- (i-1)*n_samples + seq_len(n_samples)
      idx01              <- (N_bits + i-1)*n_samples + seq_len(n_samples)
      fr[idx10, "value"] <- qc$fliprates[,paste0("rate10_bit", i) == fr_names]
      fr[idx10, "bit"]   <- i
      fr[idx10, "type"]  <- "1>0"
      fr[idx01, "value"] <- qc$fliprates[,paste0("rate01_bit", i) == fr_names]
      fr[idx01, "bit"]   <- i
      fr[idx01, "type"]  <- "0>1"
    }
    fr$bit <- as.factor(fr$bit)
    # Make plot
    plt <- ggplot(fr, aes(bit, value, fill = type)) +
      geom_point() +
      labs(title = "Estimated Flip Rates", x = "Bit", y = "Flip Rate", color = "Flip Type") +
      theme_minimal() +
      facet_grid(type ~ .)
    return(plt)
  }
