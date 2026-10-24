
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
    bit_org,
    maxHam
  ) {
    # Check bit size 
    if (is.null(dim(codebook)))                      stop("codebook must be a matrix or data frame.")
    N_bits <- ncol(codebook)
    if (N_bits > 64)                                 stop("Codebook has more than the max-allowed 64 bits.")
    if (!is.null(bit_org)) {
      if (is.null(dim(bit_org)))                     stop("bit_org must be a matrix or data frame.")
      if (nrow(bit_org) != N_bits)                   stop("Bit organization provided, but does not match the number of bits in the codebook.")
    }
   
    # Prep codebook and STdata
    if (is.null(dim(STdata)))                        stop("STdata must be a matrix or data frame.") 
    # ... get species names
    species_names <- rownames(STdata)
    if (is.null(species_names))                      stop("STdata must have row names as barcode species names")
    # ... align rows of codebook to STdata
    if (is.null(rownames(codebook)))                 stop("codebook must have row names as barcode species names")
    if (!all(species_names %in% rownames(codebook))) stop("All barcode species in STdata must be present in codebook")
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
    
    # Check STdata columns
    req_cols <- c("rates", "variance", "counts")
    if (is.null(colnames(STdata)))                   stop("STdata must have named columns.")
    if (!all(req_cols %in% colnames(STdata)))        stop("STdata must have columns 'rates', 'variance', and 'counts'.")
   
    # Check and set maxHam (max correctable Hamming distance)
    codebook_distances <- unique_Hamming_cb(as.matrix(codebook))
    if (is.null(maxHam)) {
      maxHam <- min(codebook_distances) - 1
    } else if (!is.numeric(maxHam) || length(maxHam) != 1) {
      stop("maxHam must be a numeric of length 1.")
    } else if (maxHam >= min(codebook_distances)) {
      stop(paste0("maxHam must be less than the minimum Hamming distance between codebook entries (", min(codebook_distances), ")"))
    }
    
    return(
      list(
        STdata   = STdata,
        codebook = codebook,
        bit_org  = bit_org,
        maxHam   = maxHam
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

#' Helper function to build a tiered Gaussian prior on bit-flip correlations from imaging-round plus channel structure
#'
#' @inheritParams misread.qc
#' @param fliprate_priors A \code{fliprate_priors} list (see the "Flip-rate priors" section below for the full element list). If it does not contain \code{sd_round}, \code{sd_channel}, or \code{sd_other}, default values (of 0.2, 0.1, and 0.01, respectively) are provided.
#' @return The original input \code{fliprate_priors}, minus \code{sd_round}, \code{sd_channel}, and \code{sd_other} (if provided), and enriched with \code{latent_bit_round}, \code{latent_bit_channel}, and \code{sdcorr_pairs} (if \code{bit_org} supplied).
#' @inheritSection misread.qc Flip-rate priors
#' @export
build_corr_priors <- function(
    bit_org         = NULL,
    fliprate_priors = list()
  ) {
   
    # Initialize list for flip rate priors
    fr_priors <- list()
    if (length(fliprate_priors) > 0) {
      if (class(fliprate_priors) != "list") stop("fliprate_priors must be a list")
      if (all(c("w_round", "w_channel") %in% names(fliprate_priors))) {
        if (fliprate_priors$w_round + fliprate_priors$w_channel > 1) stop("w_round + w_channel must be <= 1")
      }
      fr_priors <- c(fr_priors, as.list(fliprate_priors))
    }
    N_bits <- 0
    
    # Extract bit structure, integer-coded (via factor()) so arbitrary label types (character, numeric, ...) work
    if (!is.null(bit_org)) {
      if (is.null(dim(bit_org)))                                   stop("bit_org must be a matrix or data frame with N_bit rows")
      if (!all(c("imagingRound", "color") %in% colnames(bit_org))) stop("bit_org must have columns imagingRound and color")
      N_bits                       <- nrow(bit_org)
      fr_priors$latent_bit_round   <- as.integer(factor(bit_org$imagingRound))
      fr_priors$latent_bit_channel <- as.integer(factor(bit_org$color))
    } else {
      return(fr_priors)
    }
   
    # Compute number of free correlation parameters
    corr_free    <- N_bits * (N_bits - 1) / 2
    sdcorr_pairs <- numeric(corr_free)
    
    # Build sdcorr_pairs
    sd_round   <- 0.2
    sd_channel <- 0.1
    sd_other   <- 0.01
    if ("sd_round" %in% names(fr_priors)) {
      sd_round             <- fr_priors$sd_round
      fr_priors$sd_round   <- NULL
    }
    if ("sd_channel" %in% names(fr_priors)) {
      sd_channel           <- fr_priors$sd_channel
      fr_priors$sd_channel <- NULL
    }
    if ("sd_other" %in% names(fr_priors)) {
      sd_other             <- fr_priors$sd_other
      fr_priors$sd_other   <- NULL
    }
    # Fill in the same order as main.cpp's strict lower triangle (pair (i0, j0), 0-based, j0 < i0,
    # at linear index i0*(i0-1)/2 + j0): looping i0 = 1..N_bits-1, j0 = 0..i0-1 visits that index in
    # increasing order, so a plain running counter k matches it without computing the formula.
    k <- 0
    for (i0 in 1:(N_bits - 1)) {
      for (j0 in 0:(i0 - 1)) {
        k <- k  + 1
        i <- i0 + 1 # 1-based R index into bit_round/bit_channel
        j <- j0 + 1
        if (fr_priors$latent_bit_round[i] == fr_priors$latent_bit_round[j]) {
          sdcorr_pairs[k] <- sd_round
        } else if (fr_priors$latent_bit_channel[i] == fr_priors$latent_bit_channel[j]) {
          sdcorr_pairs[k] <- sd_channel
        } else {
          sdcorr_pairs[k] <- sd_other
        }
      }
    }
    fr_priors$sdcorr_pairs <- sdcorr_pairs
   
    # Return parameters
    return(fr_priors)
  }

#' Run misread QC analysis on summary stats data from spatial transcriptomics experiment
#'
#' This function takes summary statistics (\code{STdata}) from a FISH-based spatial transcriptomics experiment and the barcode codebook (including blanks labelled with "Blank") and runs L-BFGS (via nlopt) to make a best-fit estimate of the bit-flip rates and correlations while using the delta method to estimate 95% confidence intervals. The estimation uses an analytic conditional probability model to compute, for each barcode \emph{b}, values expected based on bit-flip rates and bit-flip correlations, for: \itemize{
#'   \item \strong{Read Count}: The number of spots labelled \emph{b} before error correction.
#'   \item \strong{Corrected Count}: The number of spots labelled \emph{b} after error correction.
#'   \item \strong{Hit Count}: The number of spots labelled \emph{b} after error correction which are in fact \emph{b}, i.e., the number of reads error-corrected to \emph{b} which are correct, aka a "hit". 
#' }
#'  The "best fit" is defined as least mean squared log error (msle) computed by comparing the analytically implied expected corrected counts to the observed corrected counts included in \code{STdata}. If \code{counts_uncorrected} are included in \code{STdata} and \code{weights$erc} > 0, then the msle is computed for the read count as well and combined with the corrected-count msle, weighted by \code{weights$erc}. Counts (corrected and read) for blank barcodes can be upweighted via \code{weights$blank}. Penalty terms based on flip rate priors (expected value for flip rates, expected values and dispersion for bit-flip correlations) and bit-flip correlation feasibility can be added by setting non-zero values for \code{weights$prior}, \code{weights$dispersion}, and \code{weights$feasibility}. For each barcode, the expected "confidence ratio" (CR) and expected positive predictive value (PPV) are computed, from these implied count values. The CR is a quality control metric proposed by the original developers of MERFISH (e.g., see DOI 10.1016/bs.mie.2016.03.020) defined per barcode as the read count over the corrected count, while PPV is defined as the hit count over the corrected count. That is, PPV tells us, for each barcode \emph{b}, the expected percentage of spots decoded as \emph{b} which are in fact \emph{b}. 
#'
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts", and optionally "counts_uncorrected". The "counts" column must provide the corrected count, the "counts_uncorrected" column can be used optionally to provide the read count. The column "rates" is assumed to hold the per-cell rates for counts, empirically estimated as corrected count over the number of cells. The "variance" column is assumed to hold the variance of the corrected count across cells. The matrix must have barcode names (e.g., gene or protein species) as row names, with blank species having "Blank" in their name.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}.
#' @param bit_org Optional data frame (or matrix) with one row per bit (columns \code{imagingRound} and \code{color}) giving the imaging-round plus color-channel structure of the codebook, used to enrich \code{fliprate_priors} with tiered correlations (\code{sdcorr_pairs}) via \code{build_corr_priors()}. Default \code{NULL} leaves \code{sdcorr_pairs} empty (no correlation structure). Tiered values for use with \code{bit_org} structure is set via entries \code{sd_round}, \code{sd_channel}, and \code{sd_other} in \code{fliprate_priors} (see the "Flip-rate priors" section below).
#' @param n_forks Number of process forks to use for parallelizable computation (expected-count computation of \code{misread.qc()}, spot simulation in \code{misread.qc.DG()}), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param max_flips When analytically computing expected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between 6 and 10 is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing.
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
#' @param n_restarts Number of L-BFGS multi-starts to run. Every restart's initial point is a random draw from \code{fliprate_priors}: the first restart always uses the same fixed seed (so a single run, the default, is reproducible across calls), while later restarts each use a different (but still deterministic, not system-random) seed. The best (lowest msle) result found across all restarts is kept. Default 1 (single, fixed-seed run).
#' @param maxHam Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries. This default value may or may not be too aggressive. There is a trade-off between missing spots of barcode \emph{b} misread as not \emph{b} (low \code{maxHam}) and correcting spots not of barcode \emph{b} as \emph{b} (high \code{maxHam}). Also, lower \code{maxHam} will run faster and use less memory than higher \code{maxHam}.
#' @param fliprate_priors Structured prior list on flip rates and bit-flip correlations (mirrors the C++ \code{FlipRatePriors} struct, plus possibly additional terms to construct \code{sdcorr_pairs} from \code{bit_org}); see the "Flip-rate priors" section below for the full element list.
#' @param weights List with elements \code{blank}, \code{prior}, \code{dispersion}, \code{erc}, and \code{feasibility} (defaults 1.0, 0.0, 0.0, 0.0, 0.0): \code{blank} up-weights blank barcodes' counts in the fit objective; \code{prior} and \code{dispersion} scale the Gaussian-prior and dispersion penalty (respectively) on flip rates and bit-flip correlations from \code{fliprate_priors} (0 disables); \code{erc} weights the read-count (i.e., uncorrected) msle term combined with the corrected-count msle when \code{STdata} includes \code{counts_uncorrected} (0 disables, so only the corrected-count msle is used); \code{feasibility} penalizes bit-flip correlations for requiring a large correction before they could be realized as an actual (PSD) per-barcode correlation structure (0 disables). Any element omitted from the supplied list falls back to its default.
#' @param max_fs_bc Cap on the number of (deterministically subsampled) barcodes used to estimate the feasibility penalty each time it's evaluated (per L-BFGS objective evaluation under this method; per MCMCSA step under the DG method, \code{misread.qc.DG()}), for speed on large codebooks. Ignored when \code{weights$feasibility} is 0. Default 40.
#' @return A list giving:\itemize{
#'    \item \code{STdata}: a dataframe giving the summary data from the ST run used in the estimation.
#'    \item \code{fliprates}: a labeled vector giving the estimated (i.e., best fit) flip rates and bit-flip correlations from the optimization.
#'    \item \code{erctc_plus}: a data frame giving, for each barcode, the read, corrected, and hit counts as well as \code{CR} and \code{PPV} values implied by the analytic conditional probability model, given the values in \code{fliprates}, plus a 95\% confidence interval on \code{PPV} (\code{PPV_lower} and \code{PPV_upper}) from a delta-method (Laplace) approximation -- this propagates uncertainty in the fitted flip-rate and correlation parameters through to PPV, via a Gauss-Newton information matrix truncated to its well-conditioned directions, rather than treating PPV as a simple binomial proportion. Column names are: "erc", "ecc", "ehc", "CR", "PPV", "PPV_lower", "PPV_upper".
#'    \item \code{msle}: A numeric value giving the minimal mean squared log error found by the L-BFGS optimization of the flip rates and bit-flip correlations.
#' }
#' @section Flip-rate priors:
#' \code{fliprate_priors} (and, where applicable, \code{stip_fliprate_priors}) is a structured prior
#' list on flip rates and bit-flip correlations (mirrors the C++ \code{FlipRatePriors} struct, plus
#' possibly additional terms to construct \code{sdcorr_pairs} from \code{bit_org}), with
#' elements:\itemize{
#'    \item \code{expected10}: Expected rate at which 1 bits flip to 0s, default 0.01.
#'    \item \code{expected01}: Expected rate at which 0 bits flip to 1s, default 0.05.
#'    \item \code{sd10}: Standard deviation (sd) of the 1-to-0 flip rate, default 0.1.
#'    \item \code{sd01}: sd of the 0-to-1 flip rate, default 0.2.
#'    \item \code{expectedcorr}: Expected baseline flip-correlation between bits, before considering channels and rounds, default 0.0.
#'    \item \code{sdcorr}: sd of flip-correlation between bits, before considering channels and rounds, default 0.2.
#'    \item \code{round_corr_sd}: sd fed into each round's own local correlation draw, default 0.2.
#'    \item \code{channel_rho}: Equicorrelation value shared by every pair within a channel, before w_channel scaling (must be >= -1/(m-1) for the smallest channel group of size m), default 0.5.
#'    \item \code{w_round}: Fraction of each bit's variance attributed to round membership (0 disables the round tier), default 0.7.
#'    \item \code{w_channel}: Fraction attributed to channel membership (0 disables the channel tier); w_round + w_channel must be <= 1, default 0.1.
#'    \item \code{sd_round}: For use with bit_org; Prior sd for bit pairs sharing an imaging round (tier 1, widest/most permissive), default 0.2 (matches the package's prior \code{sdcorr} default).
#'    \item \code{sd_channel}: For use with bit_org; Prior sd for bit pairs sharing a color channel but not a round (tier 2), default 0.1.
#'    \item \code{sd_other}: For use with bit_org; Prior sd for all remaining bit pairs (tier 3, pinned near \code{expectedcorr}), default 0.01.
#' } Any element omitted from the supplied list falls back to its default. These values are used as
#' the Gaussian prior-penalty term (when \code{weights$prior} or \code{weights$dispersion} are
#' non-zero) and to draw the initial parameter point(s) before fitting begins. When \code{bit_org} is
#' also supplied, \code{latent_bit_round} and \code{latent_bit_channel} will be extracted
#' automatically and \code{sdcorr_pairs} will be automatically built from these two factors plus the
#' values for \code{sd_round}, \code{sd_channel}, and \code{sd_other} (see \code{build_corr_priors()}).
#' @export
misread.qc <- function(
    STdata,
    codebook,
    bit_org         = NULL, 
    n_forks         = 1,
    max_flips       = 0,
    report_freq     = 10,
    maxeval         = 500, 
    n_restarts      = 1,
    maxHam          = NULL,
    max_fs_bc       = 40,
    fliprate_priors = list(),
    weights         = list()
  ) {
    cat("\nRunning misread QC with L-BFGS (nlopt)")
    cat("\nMax evaluations:", maxeval)
    
    # Prep inputs
    n_forks         <- check_forks(n_forks)                                   # Confirm forking settings
    data_check      <- check_data(STdata, codebook, bit_org, maxHam)          # Prep data
    fliprate_priors <- build_corr_priors(data_check$bit_org, fliprate_priors) # Build flip rate priors from bit_org
   
    # Run misread QC algorithm with L-BFGS (or MCMCSA)
    qc <- mQC(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(data_check$maxHam),
      as.integer(n_forks),
      as.integer(max_flips),
      as.integer(report_freq),
      as.integer(maxeval),
      as.list(fliprate_priors),
      as.list(weights),
      as.integer(n_restarts),
      as.integer(max_fs_bc)
    )
    
    # Annotate result output
    cat("\nBuilding summary tables")
    N_bits                   <- ncol(codebook)
    qc$erctc_plus            <- as.data.frame(qc$erctc_plus)
    row.names(qc$erctc_plus) <- qc$STdata$species
    names(qc$fliprates)      <- c(
      paste0("rate10_bit", seq_len(N_bits)),
      paste0("rate01_bit", seq_len(N_bits)),
      paste0("corr_",      seq_len(length(qc$fliprates) - 2*N_bits))
    )
    
    return(qc)
   
  }

#' Run misread QC analysis using the dichotomized-Gaussian (DG) spot-simulation plus MCMCSA method
#'
#' \code{misread.qc.DG()} implements literal correlated-Gaussian-noise spot
#' simulations, decoded and compared against \code{STdata}'s observed barcode counts via a Markov
#' Chain Monte Carlo with Simulated Annealing (MCMCSA) search over bit-flip rates and bit-flip
#' correlations. The random-walk and annealing schedules use a two-phase FR_search-then-burn-in design, 
#' stochastic bit/pair pruning, and burn-in mean-rate enforcement.
#'
#' @inheritParams misread.qc
#' @param resamples Number of steps to take in the post-burn-in random walk (these are the steps summarized in the returned results).
#' @param burnin Number of burn-in steps before resampling begins, default \code{resamples/2}.
#' @param FR_search Number of initial steps used to search for a reasonable overall (mean) flip rate before burn-in, default \code{resamples/2}.
#' @param max_step_size Size of the largest per-step Gaussian proposal in the random walk, default 0.2.
#' @param min_step_size_divider Step size is scaled down linearly to \code{max_step_size/min_step_size_divider} over burn-in, default 10.
#' @param max_temperature Simulated-annealing temperature at the start of the walk, default 0.2.
#' @param min_temperature_scale Temperature is scaled down linearly to \code{max_temperature/(min_step_size_divider/min_temperature_scale)} over burn-in, default 0.5.
#' @param n_tracker_updates Number of progress updates to print over the full walk, default 100.
#' @param N_cells Number of cells to simulate at each step, default 1000. Lowering this trades simulation accuracy for speed -- useful for quick/development runs. Caution: too small a value (a few dozen, on a codebook with many bits) can make the MCMCSA walk itself pathologically slow rather than just less accurate -- the per-step simulation becomes so noisy that the acceptance rate can collapse to near zero once the annealing temperature cools, requiring many thousands of proposals to advance. If a run seems stuck (no progress message for a long time), try raising \code{N_cells} before assuming something is wrong.
#' @param weights List with elements \code{max_blank}, \code{min_blank}, \code{prior}, \code{dispersion}, \code{feasibility}: \code{max_blank} and \code{min_blank} specify the weight given to blanks (vs. genes) in the fit objective, scheduled from \code{max_blank} at the start of \code{FR_search} down to \code{min_blank} (its final, held value) by the end of burn-in -- defaults are \code{N_genes/N_blanks + 1} and 2, respectively, as \code{FR_search} needs a strong blank signal early on to find a reasonable flip rate; The weights \code{prior}, \code{dispersion}, and \code{feasibility} all default to 0. The weights \code{prior} and \code{dispersion} are the same Gaussian-prior and dispersion penalty \code{misread.qc} uses (added to the simulated-annealing energy); \code{feasibility} penalizes the matrices \code{corr1} and \code{corr2} for requiring a large correction before they could be realized as an actual (PSD) per-barcode Gaussian correlation structure, mirroring \code{misread.qc}'s own \code{feasibility} weight so the two methods' fitted correlations stay comparable.
#' @return A list giving:\itemize{
#'    \item \code{sim_summaries}: a data frame, one row per accepted MCMCSA step, of simulation-fit summary statistics (CR/PPV slope and MSE, stipulated flip rates/correlations, decoding rates, etc.).
#'    \item \code{PPV}: a data frame, one row per accepted MCMCSA step and one column per gene barcode, of the simulated PPV for that barcode at that step.
#'    \item \code{PPV_summary}: a data frame summarizing \code{PPV} (over the \code{resamples} rows following burn-in) per gene barcode: median, lower/upper 95\% CI, and whether the lower CI clears 0.8.
#'    \item \code{Counts_sorted_sim} and \code{Counts_sorted_obs}: data frames, one row per step, of the simulated and observed barcode counts (sorted consistently for comparison).
#'    \item \code{PPV_plot}: a ggplot object visualizing \code{PPV_summary}.
#'    \item \code{sorted_count_plot}: a ggplot object comparing simulated vs. observed sorted barcode counts over the \code{resamples} rows following burn-in.
#' }
#' @inheritSection misread.qc Flip-rate priors
#' @export
misread.qc.DG <- function(
    STdata,
    codebook,
    bit_org                  = NULL,
    n_forks                  = 1,
    maxHam                   = NULL,
    max_fs_bc                = 40,
    resamples                = 1000,
    burnin                   = NULL,
    FR_search                = NULL,
    max_step_size            = 0.2,
    min_step_size_divider    = 10,
    max_temperature          = 0.2,
    min_temperature_scale    = 0.5,
    n_tracker_updates        = 100,
    N_cells                  = 1e3,
    fliprate_priors          = list(),
    weights                  = list()
  ) {
    cat("\n\nEstimating FPQC metrics with Markov chain Monte Carlo simulated annealing (DG method)")
   
    # Prep inputs
    n_forks         <- check_forks(n_forks)                                   # Confirm forking settings
    data_check      <- check_data(STdata, codebook, bit_org, maxHam)          # Check data
    fliprate_priors <- build_corr_priors(data_check$bit_org, fliprate_priors) # Build flip rate priors from bit_org
   
    # Resolve NULL schedule defaults (same logic misread.qc.DG.legacy() uses)
    blank_mask <- grepl("Blank", rownames(data_check$STdata), ignore.case = FALSE)
    N_genes    <- sum(!blank_mask)
    N_blanks   <- sum(blank_mask)
    if (is.null(FR_search)) FR_search <- as.integer(round(resamples / 2))
    if (is.null(burnin))    burnin    <- as.integer(round(resamples / 2))
    max_blank <- N_genes / N_blanks + 1
    if (is.null(weights$max_blank)) {
      weights$max_blank <- max_blank
    } else if (weights$max_blank > max_blank) {
      cat("\n\nRequested max blank weight too high, resetting to limit:", max_blank)
      weights$max_blank <- max_blank
    }
    if (is.null(weights$min_blank)) {
      weights$min_blank <- 2
    } else if (weights$min_blank < 1) {
      cat("\n\nRequested min blank weight too low, resetting to limit: 1")
      weights$min_blank <- 1
    }
   
    # Run the MCMCSA walk
    res <- mQC_DG(
      as.matrix(data_check$STdata),
      as.matrix(data_check$codebook),
      as.integer(data_check$maxHam),
      as.integer(resamples),
      as.integer(burnin),
      as.integer(FR_search),
      as.double(max_step_size),
      as.integer(min_step_size_divider),
      as.double(max_temperature),
      as.double(min_temperature_scale),
      as.integer(n_tracker_updates),
      as.integer(n_forks),
      as.integer(N_cells),
      as.list(fliprate_priors),
      as.list(weights),
      as.integer(max_fs_bc)
    )
   
    # Annotate results. Plotting stays in R, reusing dg_legacy.R's private .DG_plot_* helpers (same
    # package namespace; that file itself is untouched).
    gene_names               <- rownames(data_check$STdata)[res$gene_idx + 1]
    sim_summaries            <- as.data.frame(res$sim_summaries)
    colnames(sim_summaries)  <- c(
      "CR_PPV_slope", "CR_PPV_mse", "sim_mse",
      "PPV_mean_filtered", "PPV_sd_filtered", "PPV_mean_retained", "ratio_genes_CR_filtered",
      "flip_rate_10_stip", "flip_rate_10", "flip_rate_01_stip", "flip_rate_01",
      "bit_noise_cor_bias_stip", "bit_noise_cor_sd_stip",
      "spot_decoding_rate", "spot_decoding_rate_corrected"
    )
    PPV               <- as.data.frame(res$PPV)
    colnames(PPV)     <- gene_names
    Counts_sorted_sim <- as.data.frame(res$Counts_sorted_sim)
    Counts_sorted_obs <- as.data.frame(res$Counts_sorted_obs)
   
    return(list(
      sim_summaries     = sim_summaries,
      PPV               = PPV,
      PPV_summary       = PPV_summary,
      Counts_sorted_sim = Counts_sorted_sim,
      Counts_sorted_obs = Counts_sorted_obs,
      PPV_plot          = PPV_plot,
      sorted_count_plot = sorted_count_plot
    ))
   
  }

#' Benchmark \code{misread.qc} function with Bernoulli simulations
#' 
#' This function takes the same summary statistics (\code{STdata}) and barcode codebook (\code{codebook}) as \code{misread.qc} and runs Bernoulli simulations with stipulated bit-flip rates and bit-flip correlations in order to estimate how well the L-BFGS algorithm recovers the bit-flip rates and bit-flip correlations for the given data set and codebook. 
#' 
#' @inheritParams misread.qc
#' @param n_sims Number of Bernoulli simulations to run. The default is 100.
#' @param stip_fliprate_priors Priors used ONLY for drawing the simulated "true" flip rates/correlations (the \code{fr_stip} returned value), as opposed to \code{fliprate_priors}, which is used for fitting as usual. Default \code{list()} falls back to \code{fliprate_priors} itself, i.e. truth and fit share one prior (prior package behavior). Pass a different list here -- e.g. from \code{build_corr_priors()} with a near-zero \code{sd_other}, to simulate a correlation structure (such as bit pairs sharing an imaging round/channel) -- while fitting with a different prior (e.g. an unstructured one, or one with different tier weights) via \code{fliprate_priors}. This is what lets \code{sim.benchmark} test whether a structured fitting prior actually helps when the simulated truth really has that structure, as opposed to when it doesn't.
#' @return A list giving:\itemize{
#'    \item \code{fr__est}: A matrix giving the flip rates and bit-flip correlations estimated by \code{misread.qc} (columns, named "rate10_bit*", "rate01_bit*", and "corr_*"), for each set of observed counts generated by Bernoulli simulation (rows).
#'    \item \code{fr_stip}: A vector giving the stipulated flip rates and bit-flip correlations used to generate the Bernoulli simulations. All simulations use the same stipulated values.
#'    \item \code{PPV_est}: A matrix giving positive predictive value (PPV) expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows).
#'    \item \code{PPV_expected}: A vector giving the PPV values expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}.
#'    \item \code{ecc_est}: A matrix giving the corrected count expected based on analytical computation, using the values in \code{fr_est}, per barcode (columns), for each simulation run (rows).
#'    \item \code{ecc_expected}: A vector giving the corrected count expected based on analytical computation for each barcode, given the stipulated flip rates and stipulated bit-flip correlations in \code{fr_stipulated}.
#'    \item \code{sim_counts}: A matrix giving the actual simulated count per barcode (columns) for each simulation (rows).
#'    }
#' @inheritSection misread.qc Flip-rate priors
sim.benchmark <- function(
    STdata,
    codebook,
    bit_org              = NULL,
    n_sims               = 5,
    n_forks              = 1,
    max_flips            = 0,
    report_freq          = 10,
    maxeval              = 500,
    n_restarts           = 1,
    maxHam               = NULL,
    fliprate_priors      = list(),
    stip_fliprate_priors = list(),
    weights              = list()
  ) {
    cat("\nBenchmarking misread QC with Bernoulli simulations, using L-BFGS (nlopt)")

    # Confirm forking is possible and check number of cores
    n_forks    <- check_forks(n_forks)
    
    # Check bit size, prep codebook and STdata, and set maxHam
    data_check <- check_data(STdata, codebook, bit_org, maxHam)
    
    # Build flip rate priors from bit_org
    N_bits               <- ncol(codebook)
    fliprate_priors      <- build_corr_priors(bit_org = data_check$bit_org, fliprate_priors = fliprate_priors)
    stip_fliprate_priors <- build_corr_priors(bit_org = data_check$bit_org, fliprate_priors = stip_fliprate_priors)
    if (!is.null(stip_fliprate_priors$latent_bit_round) && !is.null(stip_fliprate_priors$latent_bit_channel)) {
      stip_fliprate_priors$use_latent_structure <- TRUE
    }
   
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
      as.list(fliprate_priors),
      as.list(stip_fliprate_priors),
      as.list(weights),
      as.integer(n_restarts)
    )
   
    # Name flip-rate/correlation parameters (columns of fr__est, entries of fr_stip)
    fr_param_names           <- c(
      paste0("rate10_bit", seq_len(N_bits)),
      paste0("rate01_bit", seq_len(N_bits)),
      paste0("corr_",      seq_len(ncol(resids$fr__est) - 2*N_bits))
    )
    colnames(resids$fr__est) <- fr_param_names
    names(resids$fr_stip)    <- fr_param_names
    
    return(resids)
    
  }

#' Plot estimated PPV by barcode from \code{misread.qc} results
#'
#' This function takes the results from the misread.qc function and makes a plot of the estimated positive predictive value (PPV) for each barcode, along with its 95\% confidence interval (\code{PPV_lower}/\code{PPV_upper} from \code{erctc_plus}) shown as error bars. Barcodes are sorted by PPV value, in decreasing order, and a dashed red line indicates the minimum PPV cutoff for "good" barcodes. Barcodes above the cutoff are colored blue, while those that are not are colored black.
#'
#' @name plot.PPV
#' @rdname plot-PPV
#' @usage plot.PPV(
#'  qc,
#'  min_PPV = 0.8
#' )
#' @param qc List, results from misread.qc function.
#' @param min_PPV Numeric, PPV cutoff for "good" barcodes, defaults to 0.8.
#' @return ggplot object showing estimated PPV, with its 95\% confidence interval, for each barcode with cutoff.
#' @export
plot.PPV <- function(
    qc,
    min_PPV = 0.8
  ) {
   
    # Grab data
    PPV       <- qc$erctc_plus$PPV
    PPV_lower <- qc$erctc_plus$PPV_lower
    PPV_upper <- qc$erctc_plus$PPV_upper
   
    # Sort by decreasing values
    above_cutoff <- rep("Bad", length(PPV))
    above_cutoff[PPV_lower > min_PPV] <- "Good"
    ord <- order(PPV_lower, decreasing = TRUE)
    df <- data.frame(
      x            = seq_along(PPV),
      PPV          = PPV[ord],
      PPV_lower    = PPV_lower[ord],
      PPV_upper    = PPV_upper[ord],
      above_cutoff = above_cutoff[ord]
      )
   
    # Plot
    PPV_plot <- ggplot(df, aes(x = x, y = PPV)) +
      geom_errorbar(aes(ymin = PPV_lower, ymax = PPV_upper, color = above_cutoff), width = 0, alpha = 0.5) +
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
        title = "Estimated Precision (PPV) by Gene Barcode, with 95% CI")
   
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

# Plot simulation step information
plot.MCMCsteps <- function(
    sim_summaries,
    burnin,
    FR_search
  ) {
   
    # Add step column
    sim_summaries$step <- c(1:nrow(sim_summaries))
   
    # Set v-line points
    burnin             <- burnin + FR_search
    FR_search          <- mean(c(FR_search, FR_search + 1))
    burnin             <- mean(c(burnin, burnin + 1))
    vline_points       <- c(FR_search, burnin)
   
    # plot sim_mse
    plot_mse <- ggplot(sim_summaries, aes(x = step, y = sim_mse)) +
      geom_line() +
      geom_vline(
        xintercept = vline_points,
        color.     = "black",
        linetype   = "dashed",
        linewidth  = 0.5
      ) +
      labs(title = "MCMCSA step diagnostics", x = "step", y = "MSE (log)") +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = .DG_title_size),
        axis.title      = element_text(size = .DG_axis_size),
        axis.text.y     = element_text(size = .DG_axis_size),
        axis.text.x     = element_text(size = .DG_axis_size),
        legend.title    = element_text(size = .DG_legend_size),
        legend.text     = element_text(size = .DG_legend_size)
      )
   
    # plot flip rates
    df_flip <- rbind(
      data.frame(
        step     = sim_summaries$step,
        FlipRate = sim_summaries$flip_rate_10_stip,
        Type     = "1 to 0"),
      data.frame(
        step     = sim_summaries$step,
        FlipRate = sim_summaries$flip_rate_01_stip,
        Type     = "0 to 1")
    )
   
    # Plot
    plot_flip_rates <- ggplot(df_flip, aes(x = step, y = FlipRate, color = Type)) +
      geom_line() +
      geom_vline(
        xintercept = vline_points,
        color      = "black",
        linetype   = "dashed",
        linewidth  = 0.5
      ) +
      labs(title = "MCMCSA Flip Rates", x = "step", y = "Flip Rate") +
      theme_minimal() +
      theme(
        plot.title   = element_text(hjust = 0.5, size = .DG_title_size),
        axis.title   = element_text(size = .DG_axis_size),
        axis.text.y  = element_text(size = .DG_axis_size),
        axis.text.x  = element_text(size = .DG_axis_size),
        legend.title = element_text(size = .DG_legend_size),
        legend.text  = element_text(size = .DG_legend_size)
      )
   
    # plot bit noise correlations
    df_flip <- rbind(
      data.frame(
        step       = sim_summaries$step,
        corr_param = sim_summaries$bit_noise_cor_sd_stip,
        Type       = "sd"),
      data.frame(
        step       = sim_summaries$step,
        corr_param = sim_summaries$bit_noise_cor_bias_stip,
        Type       = "bias")
    ) 
   
    # Plot
    plot_corr_param <- ggplot(df_flip, aes(x = step, y = corr_param, color = Type)) +
      geom_line() +
      geom_vline(
        xintercept = vline_points,
        color      = "black",
        linetype   = "dashed",
        linewidth  = 0.5
      ) +
      labs(title = "MCMCSA bit-noise correlation parameters", x = "step", y = "Value") +
      theme_minimal() +
      theme(
        plot.title   = element_text(hjust = 0.5, size = .DG_title_size),
        axis.title   = element_text(size = .DG_axis_size),
        axis.text.y  = element_text(size = .DG_axis_size),
        axis.text.x  = element_text(size = .DG_axis_size),
        legend.title = element_text(size = .DG_legend_size),
        legend.text  = element_text(size = .DG_legend_size)
      )
   
    print(plot_mse)
    print(plot_flip_rates)
    print(plot_corr_param)
   
  }
