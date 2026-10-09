
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
    maxHam
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
    
    # Check and set maxHam
    codebook_distances <- unique_Hamming_cb(as.matrix(codebook))
    if (is.null(maxHam)) {
      maxHam <- min(codebook_distances) - 1
    } else if (maxHam >= min(codebook_distances)) {
      stop(paste0("maxHam must be less than the minimum Hamming distance between codebook entries (", min(codebook_distances), ")"))
    }
    
    return(
      list(
        STdata                           = STdata,
        codebook                         = codebook,
        maxHam = maxHam
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
#' @param maxHam Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries.
#' @param use_mcmcsa If TRUE, use the MCMCSA simulated-annealing optimizer instead of L-BFGS (nlopt), default FALSE.
#' @param mcmcsa_step_hi MCMCSA only: starting (largest) per-step Gaussian proposal SD, linearly decayed to \code{mcmcsa_step_lo} over \code{maxeval} steps, default 0.05.
#' @param mcmcsa_step_lo MCMCSA only: ending (smallest) per-step Gaussian proposal SD, default 0.005.
#' @param mcmcsa_temp_hi MCMCSA only: starting (largest) annealing temperature, linearly decayed to \code{mcmcsa_temp_lo} over \code{maxeval} steps, default 0.1.
#' @param mcmcsa_temp_lo MCMCSA only: ending (smallest) annealing temperature, default 0.01.
#' @param mcmcsa_seed MCMCSA only: base RNG seed (offset by restart number), default 12345.
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
    n_forks           = 1,
    max_flips         = 0,
    report_freq       = 10,
    maxeval           = 500,
    maxHam = NULL,
    fliprate_priors   = list(),
    blank_weight      = 1.0,
    prior_weight      = 0.0,
    dispersion_weight = 0.0,
    obs_erc           = numeric(0),
    erc_weight        = 0.0,
    n_restarts        = 1,
    use_mcmcsa        = FALSE,            # If TRUE, use the MCMCSA simulated-annealing optimizer instead of L-BFGS (nlopt)
    mcmcsa_step_hi    = 0.05,             # MCMCSA only: starting (largest) per-step Gaussian proposal SD, linearly decayed to mcmcsa_step_lo over maxeval steps
    mcmcsa_step_lo    = 0.005,            # MCMCSA only: ending (smallest) per-step Gaussian proposal SD
    mcmcsa_temp_hi    = 0.1,              # MCMCSA only: starting (largest) annealing temperature, linearly decayed to mcmcsa_temp_lo over maxeval steps
    mcmcsa_temp_lo    = 0.01,             # MCMCSA only: ending (smallest) annealing temperature
    mcmcsa_seed       = 12345             # MCMCSA only: base RNG seed (offset by restart number)
  ) {
    cat("\nRunning misread QC with", if (use_mcmcsa) "MCMCSA" else "L-BFGS (nlopt)")
    cat("\nMax evaluations:", maxeval)
    
    # Confirm forking is possible and check number of cores
    n_forks    <- check_forks(n_forks)
    
    # Check bit size, prep codebook and STdata, and set maxHam
    data_check <- check_data(STdata, codebook, maxHam)
    
    # Run misread QC algorithm with L-BFGS (or MCMCSA)
    qc <- mQC(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(data_check$maxHam),
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      as.integer(maxeval),
      fliprate_priors,
      as.double(blank_weight),
      as.double(prior_weight),
      as.double(dispersion_weight),
      as.double(obs_erc),
      as.double(erc_weight),
      as.integer(n_restarts),
      as.logical(use_mcmcsa),
      as.double(mcmcsa_step_hi),
      as.double(mcmcsa_step_lo),
      as.double(mcmcsa_temp_hi),
      as.double(mcmcsa_temp_lo),
      as.integer(mcmcsa_seed)
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

#' Benchmark \code{misread.qc} function with Bernoulli simulations
#' 
#' This function takes the same summary statistics (\code{STdata}) and barcode codebook (\code{codebook}) as \code{misread.qc} and runs Bernoulli simulations with stipulated bit-flip rates and bit-flip correlations in order to estimate how well the L-BFGS algorithm recovers the bit-flip rates and bit-flip correlations for the given data set and codebook. 
#' 
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}. 
#' @param n_sims Number of Bernoulli simulations to run. The default is 100. 
#' @param n_forks Number of process forks to use for expected-count computation (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param max_flips When analytically computing expected corrected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between six and ten is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing. 
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
#' @param maxHam Maximum Hamming distance for misreads to be corrected (max correctable Hamming distance), default NULL sets it to one less than the minimum Hamming distance between codebook entries.
#' @param use_mcmcsa If TRUE, use the MCMCSA simulated-annealing optimizer instead of L-BFGS (nlopt), default FALSE.
#' @param mcmcsa_step_hi MCMCSA only: starting (largest) per-step Gaussian proposal SD, linearly decayed to \code{mcmcsa_step_lo} over \code{maxeval} steps, default 0.05.
#' @param mcmcsa_step_lo MCMCSA only: ending (smallest) per-step Gaussian proposal SD, default 0.005.
#' @param mcmcsa_temp_hi MCMCSA only: starting (largest) annealing temperature, linearly decayed to \code{mcmcsa_temp_lo} over \code{maxeval} steps, default 0.1.
#' @param mcmcsa_temp_lo MCMCSA only: ending (smallest) annealing temperature, default 0.01.
#' @param mcmcsa_seed MCMCSA only: base RNG seed (offset by restart number), default 12345.
#' @return A list giving:\itemize{
#'    \item \code{fr__est}: A matrix giving the flip rates and bit-flip correlations estimated by \code{misread.qc} (columns, named "rate10_bit*", "rate01_bit*", and "corr_*"), for each set of observed counts generated by Bernoulli simulation (rows).
#'    \item \code{fr_stip}: A vector giving the stipulated flip rates and bit-flip correlations used to generate the Bernoulli simulations. All simulations use the same stipulated values.
#'    \item \code{PPV_est}: A matrix giving positive predictive value (PPV) expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows). 
#'    \item \code{PPV_expected}: A vector giving the PPV values expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}. 
#'    \item \code{ecc_est}: A matrix giving the corrected count expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows).  
#'    \item \code{ecc_expected}: A vector giving the corrected count expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}. 
#'    \item \code{sim_counts}: A matrix giving the actual simulated count per barcode (columns) for each simulation (rows). 
#'    }
sim.benchmark <- function(
    STdata,
    codebook,
    n_sims            = 100,
    n_forks           = 1,
    max_flips         = 0,
    report_freq       = 10,
    maxeval           = 500,
    maxHam = NULL,
    fliprate_priors   = list(),
    blank_weight      = 1.0,
    prior_weight      = 0.0,
    dispersion_weight = 0.0,
    erc_weight        = 0.0,
    n_restarts        = 1,
    use_mcmcsa        = FALSE,            # If TRUE, use the MCMCSA simulated-annealing optimizer instead of L-BFGS (nlopt)
    mcmcsa_step_hi    = 0.05,             # MCMCSA only: starting (largest) per-step Gaussian proposal SD, linearly decayed to mcmcsa_step_lo over maxeval steps
    mcmcsa_step_lo    = 0.005,            # MCMCSA only: ending (smallest) per-step Gaussian proposal SD
    mcmcsa_temp_hi    = 0.1,              # MCMCSA only: starting (largest) annealing temperature, linearly decayed to mcmcsa_temp_lo over maxeval steps
    mcmcsa_temp_lo    = 0.01,             # MCMCSA only: ending (smallest) annealing temperature
    mcmcsa_seed       = 12345             # MCMCSA only: base RNG seed (offset by restart number)
  ) {
    cat("\nBenchmarking misread QC with Bernoulli simulations, using", if (use_mcmcsa) "MCMCSA" else "L-BFGS (nlopt)")
    
    # Confirm forking is possible and check number of cores
    n_forks    <- check_forks(n_forks)
    
    # Check bit size, prep codebook and STdata, and set maxHam
    data_check <- check_data(STdata, codebook, maxHam)
    
    # Run misread QC algorithm with L-BFGS (or MCMCSA)
    resids <- test_fr_recovery(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(n_sims), 
      as.integer(data_check$maxHam),
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      as.integer(maxeval),
      fliprate_priors,
      as.double(blank_weight),
      as.double(prior_weight),
      as.double(dispersion_weight),
      as.double(erc_weight),
      as.integer(n_restarts),
      as.logical(use_mcmcsa),
      as.double(mcmcsa_step_hi),
      as.double(mcmcsa_step_lo),
      as.double(mcmcsa_temp_hi),
      as.double(mcmcsa_temp_lo),
      as.integer(mcmcsa_seed)
    )
    
    # Name flip-rate/correlation parameters (columns of fr__est, entries of fr_stip)
    N_bits                       <- ncol(codebook)
    fr_param_names                <- c(
      paste0("rate10_bit", seq_len(N_bits)),
      paste0("rate01_bit", seq_len(N_bits)),
      paste0("corr_",      seq_len(ncol(resids$fr__est) - 2*N_bits))
    )
    colnames(resids$fr__est)     <- fr_param_names
    names(resids$fr_stip)        <- fr_param_names
    
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

#' Plot estimated flip rates from \code{misread.qc} or \code{sim.benchmark} results
#' 
#' This function takes the results from either the \code{misread.qc} function (a single point estimate of the bit-flip rates) or the \code{sim.benchmark} function (bit-flip rates estimated across simulation replicates) and makes a plot showing the estimated bit-flip rates for each bit, separated by flip type (1>0 vs 0>1). When given \code{sim.benchmark} output, the plot shows the distribution (violin plots) of estimates across replicates; when given \code{misread.qc} output, it shows the single point estimate for each bit.
#' 
#' @name plot.fr
#' @rdname plot-fr
#' @usage plot.fr(qc)
#' @param qc List, results from either the \code{misread.qc} function or the \code{sim.benchmark} function.
#' @return ggplot object showing estimated flip rates for each bit and flip type; a distribution (violin plots) when \code{qc} comes from \code{sim.benchmark}, or point estimates when \code{qc} comes from \code{misread.qc}.
#' @export
plot.fr <- function(
    qc
  ) {
    # Identify source and extract flip-rate/correlation matrix (rows = samples, 
    # columns = parameters). sim.benchmark gives multiple estimates (one row per 
    # simulation replicate) in $fr__est; misread.qc gives a single point estimate 
    # as a named vector in $fliprates, treated here as a one-row matrix.
    if (!is.null(qc$fr__est)) {
      fr_mat   <- qc$fr__est
      fr_names <- colnames(fr_mat)
      if (is.null(fr_names)) stop("qc$fr__est must have column names identifying each flip-rate/correlation parameter (see sim.benchmark)")
    } else if (!is.null(qc$fliprates)) {
      fr_names <- names(qc$fliprates)
      if (is.null(fr_names)) stop("qc$fliprates must have names identifying each flip-rate/correlation parameter (see misread.qc)")
      fr_mat   <- matrix(qc$fliprates, nrow = 1, dimnames = list(NULL, fr_names))
    } else {
      stop("qc must be output from misread.qc (with $fliprates) or sim.benchmark (with $fr__est)")
    }
    
    # Make masks
    mask10    <- grepl("rate10", fr_names)
    mask01    <- grepl("rate01", fr_names)
    n_samples <- nrow(fr_mat)
    N_bits    <- sum(mask10)
    if (sum(mask01) != N_bits) stop("Number of rate10 and rate01 entries must be the same")
    # Grab data
    fr           <- matrix(NA, nrow = 2*n_samples*N_bits, ncol = 3)
    colnames(fr) <- c("value", "bit", "type")
    fr           <- as.data.frame(fr)
    for (i in seq_len(N_bits)) {
      # Indexing stipulated by the misread.qc / sim.benchmark functions
      idx10              <- (i-1)*n_samples + seq_len(n_samples)
      idx01              <- (N_bits + i-1)*n_samples + seq_len(n_samples)
      fr[idx10, "value"] <- fr_mat[, paste0("rate10_bit", i) == fr_names]
      fr[idx10, "bit"]   <- i
      fr[idx10, "type"]  <- "1>0"
      fr[idx01, "value"] <- fr_mat[, paste0("rate01_bit", i) == fr_names]
      fr[idx01, "bit"]   <- i
      fr[idx01, "type"]  <- "0>1"
    }
    fr$bit <- as.factor(fr$bit)
    # Make plot: violin distributions across replicates for sim.benchmark output, 
    # points for a single misread.qc point estimate. geom_point() doesn't render 
    # `fill`, so it's mapped to `color` instead in that case.
    if (n_samples > 1) {
      plt <- ggplot(fr, aes(bit, value, fill = type)) +
        geom_violin() +
        labs(title = "Estimated Flip Rate Distributions Across Simulations", x = "Bit", y = "Flip Rate", fill = "Flip Type")
    } else {
      plt <- ggplot(fr, aes(bit, value, color = type)) +
        geom_point(size = 3) +
        labs(title = "Estimated Flip Rates", x = "Bit", y = "Flip Rate", color = "Flip Type")
    }
    plt <- plt +
      theme_minimal() +
      facet_grid(type ~ .)
    return(plt)
  }

# #########################################################################################################
# #########################################################################################################
# HISTORICAL CODE, RESTORED VERBATIM FOR REFERENCE -- NOT WIRED UP, NOT ACTIVE.
#
# Everything between here and the matching closing brace below is a literal copy-paste of the R-level
# wrappers for the DG-model and MCMCSA-optimizer code that existed in this package's git history before
# being removed, restored on request so it doesn't have to be dug back out of `git log`/`git show`. It is
# wrapped in `if (FALSE) { ... }` so the package keeps working as-is; none of it is currently active, and
# roxygen '#'' doc markers below have been flattened to plain '#' comments so `devtools::document()`
# doesn't try to export/document these (their C++ counterparts are themselves disabled -- see src/main.cpp).
#
# Important: these two blocks were NEVER in the repo at the same time as each other, and NEVER wired
# together -- see the longer note in src/main.cpp for the full history. In short: misreadQC() (below) is
# the last-intact R wrapper (commit f1a1a6c) for MCMCSA as an alternative *optimizer* fitting the
# *analytic* model (calls the old mQC() C++ export, itself calling MCMCSA() internally) -- no DG
# simulation involved. dichot.guass.benchmark() (below) is the last-intact R wrapper (commit f4e303e) for
# the DG spot-simulator, used only as a benchmarking tool (via test_fr_recovery()) to check how well
# L-BFGS recovers stipulated parameters from DG-simulated data -- MCMCSA was long gone by that point.
#
# Both call C++ exports (mQC(), test_fr_recovery()) by their old signatures, which do not match the
# current exports of the same/similar names in src/main.cpp -- expect to reconcile this by hand.
# #########################################################################################################
# #########################################################################################################
if (FALSE) {

# ===== restored from commit f1a1a6c: misreadQC() (R wrapper for the MCMCSA optimizer) =====

# This function takes summary statistics from a FISH-based spatial transcriptomics experiment and the barcode codebook (including blanks labelled with "Blank") and runs a Markov Chain Monte Carlo with Coupled Simulated Annealing (MCMCSA) algorithm to estimate the bit-flip rates and correlations, as well as the expected read, error-corrected, and true counts for each barcode. The function returns a list of these estimates across all iterations of the MCMCSA walk, as well as summary statistics on the flip rates and error-corrected counts. The aim is to compute both the "confidence ratio" (CR) and positive predictive value (PPV) for each barcode. 
#
# @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts", must have barcode names as row names
# @param codebook Codebook with row names as barcodes and columns as bits, must have barcode names as row names
# @param max_fr Maximum flip rate to consider in the MCMCSA algorithm, default 0.1
# @param max_corr Bit-flip correlations have lower and upper bounds of -max_corr and max_corr, default is 0.2
# @param rate10_scale Assume that 1>0 flips occur in this proportion to 0>1 flips, default is 0.2
# @param initial_corr Initial max absolute value for bit-flip correlation in the MCMCSA algorithm, default is 0.01
# @param n_steps Number of steps to run the MCMCSA algorithm, default is 1000
# @param n_forks Number of parallel forks to use for MCMCSA, default is 1 (must be 1 for Windows, can be >1 for Linux/Mac)
# @param step_size_range Numeric vector of length 2, giving the max and min step size for the MCMCSA algorithm, which will be decayed linearly over n_steps, defaults to c(0.05, 0.005)
# @param temp_range Numeric vector of length 2, giving the max and min temperature for the MCMCSA algorithm, which will be decayed linearly over n_steps, defaults to c(0.1, 0.01)
# @param corr_step_scale Numeric, giving the scale of the step size for bit-flip correlations in the MCMCSA algorithm relative to the step size for flip rates, default is 0.1
# @param maxHam Maximum Hamming distance for misreads to be corrected, default is NULL which will set it to one less than the minimum Hamming distance between codebook entries
# @param ran_seed Random seed for MCMCSA algorithm, default is 12345
# @return List giving \code{STdata}, a dataframe giving the summary data from the ST run used in the simulation, \code{fliprates}, a matrix giving the estimated flip rates and bit-flip correlations from each iteration of the MCMCSA algorithm, \code{erc}, \code{ecc}, and \code{etc}, matrices giving the estimated expected read, error-corrected, and true (i.e., correctly corrected) counts for each barcode at each iteration of the MCMCSA walk, \code{CR} and \code{PPV}, matrices giving estimated confidence ratio and positive predictive values for each iteration of the MCMCSA walk, and \code{fliprates_summary} and \code{bc_summary}, which give summary statistics on the flip rates and error-corrected counts across all iterations of the MCMCSA algorithm. 
# @export
misreadQC <- function(
    STdata, 
    codebook, 
    max_fr = 0.1,
    max_corr = 0.2,
    rate10_scale = 0.2,
    initial_corr = 0.01,
    n_steps = 1000,
    n_forks = 1,
    step_size_range = c(0.05, 0.005), 
    temp_range = c(0.1, 0.01), 
    corr_step_scale = 0.1,
    maxHam = NULL,
    ran_seed = 12345
  ) {
    cat("\nRunning misread QC with MCMCSA")
    cat("\nMax flip rate:", max_fr)
    cat("\nNumber of steps:", n_steps)
    
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
    
    # Check and set maxHam
    codebook_distances <- unique_Hamming_cb(as.matrix(codebook))
    if (is.null(maxHam)) {
      maxHam <- min(codebook_distances) - 1
    } else if (maxHam >= min(codebook_distances)) {
      stop(paste0("maxHam must be less than the minimum Hamming distance between codebook entries (", min(codebook_distances), ")"))
    }
    
    # Run misread QC algorithm with MCMCSA
    qc <- mQC(
      as.matrix(STdata), 
      as.matrix(codebook), 
      maxHam,
      c(max(step_size_range), -(max(step_size_range) - min(step_size_range))/n_steps, min(step_size_range)), # step size, initial, slope, min
      c(max(temp_range), -(max(temp_range) - min(temp_range))/n_steps, min(temp_range)), # temp, initial, slope, min
      max_fr,
      max_corr,
      initial_corr,
      corr_step_scale,
      rate10_scale,
      n_steps,
      n_forks,
      ran_seed
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
    step_range <- c(round(n_steps/2):n_steps) # take second half of MCMCSA walk to compute means and CIs
    for (p in qc_names) {
      if (p == "STdata" || p == "msle") next
      if (n_steps == 1) {
        p_means <- qc[[p]]
        ci <- rbind(p_means, p_means)
      } else {
        p_means <- colMeans(qc[[p]][step_range,])
        ci <- apply(qc[[p]][step_range,], 2, quantile, probs = c(0.025, 0.975))
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
    qc[["bc_summary"]] <- bc
    
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

# ===== restored from commit f4e303e: dichot.guass.benchmark() (R wrapper for the DG spot-simulator) =====


# Benchmark \code{misread.qc} function with dichotomized-Gaussian simulations
# 
# This function takes the same summary statistics (\code{STdata}) and barcode codebook (\code{codebook}) as \code{misread.qc} and runs dichotomized-Gaussian simulations with stipulated bit-flip rates and bit-flip correlations in order to estimate how well the L-BFGS algorithm recovers the bit-flip rates and bit-flip correlations for the given data set and codebook. 
# 
# @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names.
# @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}. 
# @param n_sims Number of dichotomized-Gaussian simulations to run. The default is 100. 
# @param n_forks Number of process forks to use for expected-count computation (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
# @param max_flips When analytically computing expected corrected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between six and ten is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing. 
# @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
# @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
# @param maxHam Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries.
# @return A list giving:\itemize{
#    \item \code{fr_est}: A matrix giving the flip rates and bit-flip correlations estimated by \code{misread.qc} (columns), for each set of observed counts generated by dichotomized-Gaussian simulation (rows).
#    \item \code{fr_stipulated}: A vector giving the stipulated flip rates and bit-flip correlations used to generate the dichotomized-Gaussian simulations. All simulations use the same stipulated values.
#    \item \code{PPV_est}: A matrix giving positive predictive value (PPV) expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows). 
#    \item \code{PPV_expected}: A vector giving the PPV values expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}. 
#    \item \code{ecc_est}: A matrix giving the corrected count expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows).  
#    \item \code{ecc_expected}: A vector giving the corrected count expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}. 
#    \item \code{sim_counts}: A matrix giving the actual simulated count per barcode (columns) for each simulation (rows). 
#    }
dichot.guass.benchmark <- function(
    STdata,
    codebook,
    n_sims      = 100,
    n_forks     = 1,
    max_flips   = 0,
    report_freq = 10,
    maxeval     = 500,
    maxHam      = NULL
  ) {
    cat("\nBenchmarking misread QC with dichotomized-Gaussian simulation")
    
    # Confirm forking is possible and check number of cores
    n_forks    <- check_forks(n_forks)
    
    # Check bit size, prep codebook and STdata, and set maxHam
    data_check <- check_data(STdata, codebook, maxHam)
    
    # Run misread QC algorithm with L-BFGS
    resids <- test_fr_recovery(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(n_sims), 
      as.integer(data_check$maxHam),
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      as.integer(maxeval),
      list()
    )
    
    return(resids)
    
  }


} # end restored historical DG/MCMCSA code
