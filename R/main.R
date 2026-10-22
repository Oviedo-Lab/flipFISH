
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
    
    # Check STdata columns 
    req_cols <- c("rates", "variance", "counts")
    if (is.null(colnames(STdata)))                   stop("STdata must have named columns.") 
    if (!all(req_cols %in% colnames(STdata)))        stop("STdata must have columns 'rates', 'variance', and 'counts'.")
    
    # Split STdata
    obs_erc <- numeric(0)
    if ("counts_uncorrected" %in% colnames(STdata)) {
      obs_erc <- STdata[,"counts_uncorrected"]
      STdata  <- STdata[,"counts_uncorrected" != colnames(STdata)]
    }
    
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
        obs_erc  = obs_erc,
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

#' Build a tiered Gaussian prior on bit-flip correlations from imaging-round plus channel structure
#'
#' \code{misread.qc} and \code{sim.benchmark} regularize the (otherwise weakly identified -- see the
#' package's benchmarking notes tutorial) pairwise bit-flip correlation parameters with a single Gaussian
#' prior shared by all bit pairs (\code{fliprate_priors$expectedcorr} and \code{fliprate_priors$sdcorr}). 
#' This function instead builds a \emph{tiered} version of that same prior: bit pairs that share an imaging 
#' round get one (wider, i.e. more permissive) prior standard deviation (sd), pairs that only share a color 
#' channel get a second (narrower) sd, and all remaining pairs get a third, tight, near-zero sd. The rationale 
#' is that correlated misreads in a MERFISH/MERSCOPE-style experiment most plausibly arise from per-round
#' effects (shared hybridization/imaging cycle) or, secondarily, per-channel effects (shared optical path across
#' rounds); pairs sharing neither have no obvious mechanism to be correlated and the data cannot usefully constrain 
#' them anyway, so they are pulled tightly to \code{expectedcorr} instead of being estimated as if they were as 
#' free as the rest.
#'
#' The returned list is a complete, ready-to-use \code{fliprate_priors} argument for \code{misread.qc} and 
#' \code{sim.benchmark} (it includes the usual \code{expected10}/\code{expected01}/\code{sd10}/\code{sd01}/
#' \code{expectedcorr}/\code{sdcorr} fields alongside the new \code{sdcorr_pairs} vector, so nothing else needs 
#' to be merged in).
#'
#' The same list can also be passed as \code{sim.benchmark}'s \code{stip_fliprate_priors} argument, to
#' simulate a "true" correlation structure organized by round/channel instead of i.i.d. noise -- but
#' there, \code{sdcorr_pairs} is only an \emph{approximation} of that structure: \code{draw_random_fr}
#' still draws each pair independently and runs the result through its usual global Cholesky
#' normalization, which leaks some of the round/channel tiers' magnitude into the pairs meant to be
#' near zero. For an exact, non-leaky version of that organized structure, use \code{build_structured_corrs()} instead.
#' 
#' @param N_bits Integer, number of bits in the codebook.
#' @param bit_round Vector of length \code{N_bits} giving each bit's imaging-round label (any type that supports \code{==}, e.g. the \code{imagingRound} column of a MERSCOPE \code{codebook_organization.csv}, in bit order 1..N_bits).
#' @param bit_channel Vector of length \code{N_bits} giving each bit's color-channel label (e.g. the \code{color} column of the same file), in bit order 1..N_bits.
#' @param sd_round Prior standard deviation (sd) for bit pairs sharing an imaging round (tier 1, widest/most permissive), default 0.2 (matches the package's prior \code{sdcorr} default).
#' @param sd_channel Prior sd for bit pairs sharing a color channel but not a round (tier 2), default 0.1.
#' @param sd_other Prior sd for all remaining bit pairs (tier 3, pinned near \code{expectedcorr}), default 0.01.
#' @param expectedcorr Prior mean for all tiers, default 0.0.
#' @param sdcorr Fallback scalar sd, carried through unchanged for use anywhere \code{sdcorr_pairs} isn't consulted (e.g. \code{draw_random_fr}'s initial sampling); does not affect the tiered penalty itself. Default 0.2.
#' @param expected10 Passed straight through to the returned list, default 0.01 (package default).
#' @param expected01 Passed straight through to the returned list, default 0.05 (package default).
#' @param sd10 Passed straight through to the returned list, default 0.1 (package default).
#' @param sd01 Passed straight through to the returned list, default 0.2 (package default).
#' @return A named list suitable as \code{misread.qc}/\code{sim.benchmark}'s \code{fliprate_priors} argument, with \code{sdcorr_pairs} added (length \code{N_bits*(N_bits-1)/2}, in the same bit-pair order as \code{fr_stip}/\code{fr__est}'s correlation columns).
#' @export
build_corr_priors <- function(
    # tier values
    sd_round           = 0.2,
    sd_channel         = 0.1,
    sd_other           = 0.01,
    # bit organization 
    bit_org            = NULL,
    latent_bit_round   = NULL,
    latent_bit_channel = NULL,
    fliprate_priors    = list()
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
      if (!is.null(latent_bit_round) || !is.null(latent_bit_channel)) stop("Cannot provide both bit_org and latent_bit_round or latent_bit_channel")
      if (is.null(dim(bit_org)))                                      stop("bit_org must be a matrix or data frame with N_bit rows")
      if (!all(c("imagingRound", "color") %in% colnames(bit_org)))    stop("bit_org must have columns imagingRound and color")
      if (any(c("latent_bit_round", "latent_bit_channel") %in% names(fliprate_priors))) 
        stop("Cannot provide both latent_bit_round or latent_bit_channel in fliprate_priors and bit_org")
      N_bits             <- nrow(bit_org)
      latent_bit_round   <- as.integer(factor(bit_org$imagingRound))
      latent_bit_channel <- as.integer(factor(bit_org$color))
    } else if (!is.null(latent_bit_round) && !is.null(latent_bit_channel)) {
      if (length(latent_bit_round) != length(latent_bit_channel))     stop("latent_bit_round and latent_bit_channel must be the same length")
      if (any(c("latent_bit_round", "latent_bit_channel") %in% names(fliprate_priors))) 
        stop("Cannot provide latent_bit_round or latent_bit_channel independently and in fliprate_priors")
      N_bits             <- length(latent_bit_round)
      latent_bit_round   <- as.integer(factor(latent_bit_round))
      latent_bit_channel <- as.integer(factor(latent_bit_channel))
    } else {
      return(fr_priors)
    }
    
    # Load bit structure
    fr_priors$latent_bit_round   <- latent_bit_round
    fr_priors$latent_bit_channel <- latent_bit_channel
   
    # Compute number of free correlation parameters
    corr_free    <- N_bits * (N_bits - 1) / 2
    sdcorr_pairs <- numeric(corr_free)
   
    # Build sdcorr_pairs
    #   Fill in the same order as main.cpp's strict lower triangle (pair (i0, j0), 0-based, j0 < i0,
    #   at linear index i0*(i0-1)/2 + j0): looping i0 = 1..N_bits-1, j0 = 0..i0-1 visits that index in
    #   increasing order, so a plain running counter k matches it without computing the formula.
    k <- 0
    for (i0 in 1:(N_bits - 1)) {
      for (j0 in 0:(i0 - 1)) {
        k <- k  + 1
        i <- i0 + 1 # 1-based R index into bit_round/bit_channel
        j <- j0 + 1
        if (latent_bit_round[i] == latent_bit_round[j]) {
          sdcorr_pairs[k] <- sd_round
        } else if (latent_bit_channel[i] == latent_bit_channel[j]) {
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
#' This function takes summary statistics (\code{STdata}) from a FISH-based spatial transcriptomics experiment and the barcode codebook (including blanks labelled with "Blank") and runs L-BFGS (via nlopt) to make a best-fit estimate of the bit-flip rates and correlations. The estimation (i.e., L-BFGS optimization) uses an analytic conditional probability model to compute, for each barcode \emph{b}, values expected based on bit-flip rates and bit-flip correlations, for: \itemize{
#'   \item \strong{Read Count}: The number of spots labelled \emph{b} before error correction.
#'   \item \strong{Corrected Count}: The number of spots labelled \emph{b} after error correction.
#'   \item \strong{Hit Count}: The number of spots labelled \emph{b} after error correction which are in fact \emph{b}, i.e., the number of reads error-corrected to \emph{b} which are correct, aka a "hit". 
#' }
#'  The "best fit" is defined as least mean squared log error (msle) computed by comparing the analytically implied expected corrected counts to the observed corrected counts included in \code{STdata}. If \code{counts_uncorrected} are included in \code{STdata} and \code{weights$erc} > 0, then the msle is computed for the read count as well and combined with the corrected-count msle accounting the value set by \code{weights$erc}. Counts (corrected and read) for blank barcodes can be upweighted via \code{weights$blank}. Penalty terms based on flip rate priors (expected value for flip rates, expected values and dispersion for bit-flip correlations) and bit-flip correlation feasibility can be added by setting non-zero values for \code{weights$prior}, \code{weights$dispersion}, and \code{weights$feasibility}. For each barcode, the expected "confidence ratio" (CR) and expected positive predictive value (PPV) are computed, from these implied count values. The CR is a quality control metric proposed by the original developers of MERFISH (e.g., see DOI 10.1016/bs.mie.2016.03.020) defined per barcode as the read count over the corrected count, while PPV is defined as the hit count over the corrected count. That is, PPV tells us, for each barcode \emph{b}, the expected percentage of spots decoded as \emph{b} which are in fact \emph{b}. 
#'
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts", and optionally "counts_uncorrected". The "counts" column must provide the corrected count, the "counts_uncorrected" column can be used to optionally provide the read count. The column "rates" is assumed to hold the per-cell rates for counts, i.e., corrected count over the number of cells. The "variance" column is assumed to hold the variance of the corrected count across cells. The matrix must have barcode names (e.g., gene or protein species) as row names, with blank species having "Blank" in their name.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}.
#' @param bit_org Optional data frame/matrix with one row per bit (columns \code{imagingRound} and \code{color}) giving the imaging-round plus color-channel structure of the codebook, used to build a tiered \code{fliprate_priors} via \code{build_corr_priors()}. Default \code{NULL} skips this (unstructured prior, unless \code{fliprate_priors} is supplied directly). When \code{bit_org} is supplied, \code{fliprate_priors} may additionally include \code{build_corr_priors()}'s own \code{sd_round}/\code{sd_channel}/\code{sd_other} to tune the tiering itself (see \code{fliprate_priors} below); these are forwarded into the \code{build_corr_priors()} call, not just overlaid on its result.
#' @param n_forks Number of process forks to use for expected-count computation (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param max_flips When analytically computing expected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between six and ten is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing.
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
#' @param n_restarts Number of L-BFGS multi-starts to run. Every restart's initial point is a random draw from \code{fliprate_priors}: the first restart always uses the same fixed seed (so a single run, the default, is reproducible across calls), while later restarts each use a different (but still deterministic, not system-random) seed. The best (lowest msle) result found across all restarts is kept. Default 1 (single, fixed-seed run).
#' @param maxHam Maximum Hamming distance for misreads to be corrected, default NULL sets it to one less than the minimum Hamming distance between codebook entries. This value may or may not be too aggressive. There is a trade-off between missing spots of barcode \emph{b} misread as not \emph{b} (low \code{maxHam}) and correcting spots not of barcode \emph{b} as \emph{b} (high \code{maxHam}). Also, lower \code{maxHam} will run faster and use less memory than higher \code{maxHam}. 
#' @param fliprate_priors Structured prior list on flip rates and bit-flip correlations (mirrors the C++ \code{FlipRatePriors} struct), with elements:\itemize{
#'    \item \code{expected10}: Expected rate at which 1 bits flip to 0s, default 0.01.
#'    \item \code{expected01}: Expected rate at which 0 bits flip to 1s, default 0.05.
#'    \item \code{sd10}: Standard deviation of the 1-to-0 flip rate, default 0.1.
#'    \item \code{sd01}: Standard deviation of the 0-to-1 flip rate, default 0.2.
#'    \item \code{expectedcorr}: Expected baseline flip-correlation between bits, before considering channels and rounds, default 0.0.
#'    \item \code{sdcorr}: Standard deviation of flip-correlation between bits, before considering channels and rounds, default 0.2.
#'    \item \code{sdcorr_pairs}: Optional per-pair override of \code{sdcorr} (same strict-lower-triangle pair index as \code{corr1}/\code{corr0}, i.e. pair \code{(i,j)}, \code{j<i}, at index \code{i*(i-1)/2+j}). Empty (the default) means "not set" and falls back to \code{sdcorr}. Lets a caller express a tiered correlation prior (e.g. a wide prior for bit-pairs sharing an imaging round, a narrower one for pairs only sharing a color channel, and a tight near-zero prior for everything else) via \code{build_corr_priors()}, instead of one global sd for all pairs.
#' } Any element omitted from the supplied list falls back to its default (so \code{list()}, the function default, is equivalent to the defaults above with no tiering). These values are used as the Gaussian prior-penalty term (when \code{weights$prior}/\code{weights$dispersion} are non-zero) and to draw every restart's initial point (see \code{n_restarts}). When \code{bit_org} is also supplied, \code{fliprate_priors} may instead/additionally include \code{build_corr_priors()}'s \code{sd_round}/\code{sd_channel}/\code{sd_other} (rather than a hand-built \code{sdcorr_pairs}) to tune the per-tier sds that \code{build_corr_priors()} uses to construct \code{sdcorr_pairs} automatically; any such elements are forwarded into the \code{build_corr_priors()} call rather than overlaid on its output.
#' @param weights List with elements \code{blank}, \code{prior}, \code{dispersion}, \code{erc}, \code{feasibility} (defaults 1.0, 0.0, 0.0, 0.0, 0.0): \code{blank} upweights blank barcodes' counts in the fit objective; \code{prior}/\code{dispersion} scale the Gaussian-prior/dispersion penalty on flip rates/bit-flip correlations from \code{fliprate_priors} (0 disables); \code{erc} weights the read-count (uncorrected) msle term combined with the corrected-count msle when \code{STdata} includes \code{counts_uncorrected} (0 disables, so only the corrected-count msle is used); \code{feasibility} penalizes bit-flip correlations for requiring a large correction before they could be realized as an actual (PSD) per-barcode correlation structure (0 disables). Any element omitted from the supplied list falls back to its default.
#' @param max_feasibility_barcodes Cap on the number of (deterministically subsampled) barcodes used to estimate the feasibility penalty per objective evaluation, for speed on large codebooks. Ignored when \code{weights$feasibility} is 0. Default 40.
#' @return A list giving:\itemize{
#'    \item \code{STdata}: a dataframe giving the summary data from the ST run used in the estimation.
#'    \item \code{fliprates}: a labeled vector giving the estimated (i.e., best fit) flip rates and bit-flip correlations from the optimization. 
#'    \item \code{erctc_plus}: a data frame giving, for each barcode, the read, corrected, and hit counts as well as \code{CR} and \code{PPV} values implied by the analytic conditional probability model, given the values in \code{fliprates}, plus a 95\% confidence interval on \code{PPV} (\code{PPV_lower}/\code{PPV_upper}) from a delta-method (Laplace) approximation -- this propagates uncertainty in the fitted flip-rate/correlation parameters through to PPV, via a Gauss-Newton information matrix truncated to its well-conditioned directions, rather than treating PPV as a simple binomial proportion (a naive Wald interval on PPV alone understates its true uncertainty by roughly an order of magnitude, since it ignores parameter uncertainty entirely -- see \code{old_vs_new_ppv_comparison/} at the package root). Column names are: "erc", "ecc", "ehc", "CR", "PPV", "PPV_lower", "PPV_upper".
#'    \item \code{msle}: A numeric value giving the minimal mean squared log error found by the L-BFGS optimization of the flip rates and bit-flip correlations.
#' }
#' @export
misread.qc <- function(
    STdata,
    codebook,
    bit_org           = NULL, 
    n_forks           = 1,
    max_flips         = 0,
    report_freq       = 10,
    maxeval           = 500, 
    n_restarts        = 1,
    maxHam            = NULL,
    fliprate_priors   = list(),
    weights           = list(
      blank       = 1.0, 
      prior       = 0.0, 
      dispersion  = 0.0, 
      erc         = 0.0,
      feasibility = 0.0
    ),
    max_feasibility_barcodes = 40
  ) {
    cat("\nRunning misread QC with L-BFGS (nlopt)")
    cat("\nMax evaluations:", maxeval)
   
    # Fill in any weight the caller's list omits with its default, so a partial `weights` list
    # (e.g. missing `erc`) doesn't silently pass NULL -> numeric(0) through to mQC() below
    weights <- modifyList(
      list(blank = 1.0, prior = 0.0, dispersion = 0.0, erc = 0.0, feasibility = 0.0),
      weights
    )
   
    # Confirm forking is possible and check number of cores
    n_forks         <- check_forks(n_forks)
    # Check bit size, prep codebook and STdata, and set maxHam
    data_check      <- check_data(STdata, codebook, bit_org, maxHam)
    # Build flip rate priors from bit_org. Strip latent_bit_round/latent_bit_channel afterward: those
    # two fields switch draw_random_fr()'s initial-point draw over to the structured-truth generative
    # mode (see build_corr_priors()), which is meant for sim.benchmark()'s stip_fliprate_priors, not
    # for this function's own fitting prior -- bit_org here should only drive the tiered sdcorr_pairs
    # penalty, as before.
    N_bits                       <- ncol(codebook)
    fliprate_priors              <- build_corr_priors(bit_org = data_check$bit_org, fliprate_priors = fliprate_priors)
   
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
      as.double(weights$blank),
      as.double(weights$prior),
      as.double(weights$dispersion),
      as.double(data_check$obs_erc),
      as.double(weights$erc),
      as.integer(n_restarts),
      as.double(weights$feasibility),
      as.integer(max_feasibility_barcodes)
    )
    
    # Annotate result output
    cat("\nBuilding summary tables")
    qc$erctc_plus            <- as.data.frame(qc$erctc_plus)
    row.names(qc$erctc_plus) <- qc$STdata$species
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

#' Run misread QC analysis using the dichotomized-Gaussian (DG) spot-simulation plus MCMCSA method
#'
#' \code{misread.qc.DG} is a thin wrapper around the C++ \code{mQC_DG()} engine -- the same pattern
#' \code{misread.qc()}/\code{mQC()} use -- implementing literal correlated-Gaussian-noise spot
#' simulation, decoded and compared against \code{STdata}'s observed barcode counts via a Markov
#' Chain Monte Carlo with Simulated Annealing (MCMCSA) search over bit-flip rates and bit-flip
#' correlations. It supersedes \code{misread.qc.DG.legacy} (the verbatim-preserved historical port) in
#' three ways: (1) it correctly inverts the desired bit-flip correlations into the underlying Gaussian
#' noise correlations actually needed to simulate them -- \code{misread.qc.DG.legacy} plugged bit-flip
#' correlations in directly where Gaussian noise correlations were needed, and projects the result to
#' the nearest valid (positive-semidefinite) correlation matrix per barcode when the requested
#' correlations aren't jointly realizable; (2) it uses separate correlation matrices for when a given
#' bit reads 1 vs. 0 (\code{corr1}/\code{corr0}), matching \code{misread.qc}'s \code{FlipRates}
#' structure, instead of one shared matrix; and (3) it reuses \code{misread.qc}'s structured, tiered
#' priors (\code{fliprate_priors}/\code{bit_org}) for the MCMCSA walk's proposal step sizes and
#' acceptance-penalty term, not just its initial draw -- a tightly-pinned prior tier gets both a
#' smaller per-step proposal and a larger penalty for moving away from it. The random-walk/annealing
#' mechanics themselves (schedules, the two-phase FR_search-then-burnin design, stochastic bit/pair
#' pruning, burn-in mean-rate enforcement) are unchanged from \code{misread.qc.DG.legacy}.
#'
#' One simplification relative to \code{misread.qc.DG.legacy}: the burn-in-phase mean-flip-rate
#' enforcement targets the mean of the \emph{stipulated} (proposed/accepted) flip rates over the
#' flip-rate-search window, rather than re-deriving an independently \emph{observed} flip rate from
#' regenerating ground-truth RNA and comparing it to decoded spots. The two should track each other
#' closely in expectation; this avoids porting that secondary diagnostic machinery. Correspondingly,
#' \code{sim_summaries}' \code{flip_rate_10}/\code{flip_rate_01} columns (the observed-from-spots
#' variant) are always \code{NA} here, same as \code{misread.qc.DG.legacy} already is except near its
#' own \code{FR_search} boundary.
#'
#' @param STdata Numeric matrix or data frame with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names, with blank species having "Blank" in their name.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0. All row names from \code{STdata} must be included in the row names for \code{codebook}.
#' @param bit_org Optional data frame/matrix with one row per bit (columns \code{imagingRound}/\code{color}) giving the imaging-round/color-channel structure of the codebook, used to build a tiered \code{fliprate_priors} via \code{build_corr_priors()} -- same as \code{misread.qc}'s \code{bit_org} argument, including the ability to tune the tiering itself by also putting \code{build_corr_priors()}'s \code{sd_round}/\code{sd_channel}/\code{sd_other} in \code{fliprate_priors}. Default \code{NULL} skips this (unstructured prior, unless \code{fliprate_priors} is supplied directly).
#' @param maxHam Maximum correctable Hamming distance used for error correction, default \code{NULL} sets it to one less than the minimum Hamming distance between codebook entries.
#' @param resamples Number of steps to take in the post-burn-in random walk (these are the steps summarized in the returned results).
#' @param burnin Number of burn-in steps before resampling begins, default \code{resamples/2}.
#' @param FR_search Number of initial steps used to search for a reasonable overall (mean) flip rate before burn-in, default \code{resamples/2}.
#' @param max_step_size Size of the largest per-step Gaussian proposal in the random walk, default 0.2.
#' @param min_step_size_divider Step size is scaled down linearly to \code{max_step_size/min_step_size_divider} over burn-in, default 10.
#' @param max_temperature Simulated-annealing temperature at the start of the walk, default 0.2.
#' @param min_temperature_scale Temperature is scaled down linearly to \code{max_temperature/(min_step_size_divider/min_temperature_scale)} over burn-in, default 0.5.
#' @param n_tracker_updates Number of progress updates to print over the full walk, default 100.
#' @param n_forks Number of process forks to use for simulating spots (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param N_cells Number of cells to simulate at each step, default 1000. Lowering this trades simulation accuracy for speed -- useful for quick/development runs. Caution: too small a value (a few dozen, on a codebook with many bits) can make the MCMCSA walk itself pathologically slow rather than just less accurate -- the per-step simulation becomes so noisy that the acceptance rate can collapse to near zero once the annealing temperature cools, requiring many thousands of proposals to advance. If a run seems stuck (no progress message for a long time), try raising \code{N_cells} before assuming something is wrong.
#' @param fliprate_priors Structured prior on flip rates/bit-flip correlations, same format \code{misread.qc}/\code{build_corr_priors()} use (see \code{misread.qc}'s \code{fliprate_priors} argument for the full element list). Default \code{list()} (unstructured). Drives the MCMCSA walk's initial draw, per-pair proposal step-size tiering (via \code{sdcorr_pairs}), and the prior-penalty term (via \code{weights$prior}/\code{weights$dispersion}).
#' @param weights List with elements \code{max_blank}, \code{min_blank}, \code{prior}, \code{dispersion}, \code{feasibility}: \code{max_blank}/\code{min_blank} are the weight given to blanks (vs. genes) in the fit objective, scheduled from \code{max_blank} at the start of \code{FR_search} down to \code{min_blank} (its final, held value) by the end of burn-in -- defaults are \code{N_genes/N_blanks + 1} and 2, respectively, mirroring \code{misread.qc}'s constant \code{weights$blank} but as a schedule rather than a single scalar, since \code{FR_search} needs a strong blank signal early on to find a reasonable flip rate; \code{prior}/\code{dispersion}/\code{feasibility} all default to 0, matching \code{misread.qc}'s "0 disables it" convention: \code{prior}/\code{dispersion} are the same Gaussian-prior/dispersion penalty \code{misread.qc} uses (added to the simulated-annealing energy); \code{feasibility} penalizes corr1/corr0 for requiring a large correction before they could be realized as an actual (PSD) per-barcode Gaussian correlation structure, mirroring \code{misread.qc}'s own \code{feasibility} weight so the two methods' fitted correlations stay comparable.
#' @param max_feasibility_barcodes Cap on the number of (deterministically subsampled) barcodes used to estimate the feasibility penalty per step, for speed on large codebooks. Ignored when \code{weights$feasibility} is 0.
#' @return A list giving:\itemize{
#'    \item \code{sim_summaries}: a data frame, one row per accepted MCMCSA step, of simulation-fit summary statistics (CR/PPV slope and MSE, stipulated flip rates/correlations, decoding rates, etc.).
#'    \item \code{PPV}: a data frame, one row per accepted MCMCSA step and one column per gene barcode, of the simulated PPV for that barcode at that step.
#'    \item \code{PPV_summary}: a data frame summarizing \code{PPV} (over the \code{resamples} rows following burn-in) per gene barcode: median, lower/upper 95\% CI, and whether the lower CI clears 0.8.
#'    \item \code{Counts_sorted_sim}/\code{Counts_sorted_obs}: data frames, one row per step, of the simulated and observed barcode counts (sorted consistently for comparison).
#'    \item \code{PPV_plot}: a ggplot object visualizing \code{PPV_summary}.
#'    \item \code{sorted_count_plot}: a ggplot object comparing simulated vs. observed sorted barcode counts over the \code{resamples} rows following burn-in.
#' }
#' @export
misread.qc.DG <- function(
    STdata,
    codebook,
    bit_org                  = NULL,
    maxHam                   = NULL,
    resamples                = 1000,
    burnin                   = NULL,
    FR_search                = NULL,
    max_step_size            = 0.2,
    min_step_size_divider    = 10,
    max_temperature          = 0.2,
    min_temperature_scale    = 0.5,
    n_tracker_updates        = 100,
    n_forks                  = 1,
    N_cells                  = 1e3,
    fliprate_priors          = list(),
    weights                  = list(
      max_blank   = NULL,
      min_blank   = NULL,
      prior       = 0.0,
      dispersion  = 0.0,
      feasibility = 0.0
    ),
    max_feasibility_barcodes = 40
  ) {
   
    cat("\n\nEstimating FPQC metrics with Markov chain Monte Carlo simulated annealing (DG method)")
   
    # Fill in any weight the caller's list omits with its default (same rationale as misread.qc())
    weights <- modifyList(
      list(max_blank = NULL, min_blank = NULL, prior = 0.0, dispersion = 0.0, feasibility = 0.0),
      weights
    )
   
    # Confirm forking is possible and check number of cores (reuses misread.qc()'s own helper)
    n_forks <- check_forks(n_forks)
   
    # Check bit size, prep codebook and STdata, and set maxHam (reuses misread.qc()'s own helper;
    # strips an optional counts_uncorrected column, which mQC_DG() has no use for)
    data_check <- check_data(STdata, codebook, bit_org, maxHam)
    
    # Build flip rate priors from bit_org, same as misread.qc() -- including stripping
    # latent_bit_round/latent_bit_channel back out so bit_org only drives sdcorr_pairs here, not
    # draw_random_fr()'s structured-truth generative mode (see misread.qc()'s own comment on this)
    N_bits                       <- ncol(data_check$codebook)
    fliprate_priors              <- build_corr_priors(bit_org = data_check$bit_org, fliprate_priors = fliprate_priors)
   
    # Resolve NULL schedule defaults (same logic misread.qc.DG.legacy() uses)
    blank_mask <- grepl("Blank", rownames(data_check$STdata), ignore.case = FALSE)
    N_genes    <- sum(!blank_mask)
    N_blanks   <- sum(blank_mask)
    if (is.null(FR_search)) FR_search <- as.integer(resamples / 2)
    if (is.null(burnin))    burnin    <- as.integer(resamples / 2)
    max_blank_weight_limit <- N_genes / N_blanks + 1
    if (is.null(weights$max_blank)) {
      weights$max_blank <- max_blank_weight_limit
    } else if (weights$max_blank > max_blank_weight_limit) {
      cat("\n\nRequested max blank weight too high, resetting to limit:", max_blank_weight_limit)
      weights$max_blank <- max_blank_weight_limit
    }
    if (is.null(weights$min_blank)) {
      weights$min_blank <- 2
    } else if (weights$min_blank < 1) {
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
      as.double(weights$max_blank),
      as.double(weights$min_blank),
      as.integer(n_tracker_updates),
      as.integer(n_forks),
      as.integer(N_cells),
      as.list(fliprate_priors),
      as.double(weights$prior),
      as.double(weights$dispersion),
      as.double(weights$feasibility),
      as.integer(max_feasibility_barcodes)
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
   
    # Same resamples-row window misread.qc.DG.legacy() uses (1-indexed, inclusive of the last burn-in
    # row through the last resampled row)
    n_steps        <- nrow(PPV)
    resamples_rows <- (FR_search + burnin):(FR_search + burnin + resamples)
    resamples_rows <- resamples_rows[resamples_rows >= 1 & resamples_rows <= n_steps]
   
    PPV_plot_info <- .DG_plot_PPV_cutoff(PPV[resamples_rows, , drop = FALSE])
    PPV_summary   <- PPV_plot_info$summary_df
    PPV_plot      <- PPV_plot_info$PPV_plot
    print(PPV_plot)
   
    counts_sorted_sim_ <- apply(Counts_sorted_sim[resamples_rows, , drop = FALSE], 2, median)
    counts_sorted_obs_ <- unlist(as.vector(Counts_sorted_obs[1, ]))
    sorted_count_plot  <- .DG_plot_counts_sorted(counts_sorted_sim_, counts_sorted_obs_, length(gene_names))
    print(sorted_count_plot)
   
    .DG_plot_MCMC_steps(sim_summaries, burnin, FR_search)
   
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
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0, depending on whether an mRNA molecule of the species represented by the row is expected to luminescence in the bit represented by the column. All row names from \code{STdata} must be included in the row names for \code{codebook}. 
#' @param n_sims Number of Bernoulli simulations to run. The default is 100. 
#' @param n_forks Number of process forks to use for expected-count computation (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param max_flips When analytically computing expected corrected counts per barcode, the function will ignore misreads larger than this Hamming distance. The default is 0, which is interpreted as no limit. Using all misreads will likely be prohibitively slow; a value between six and ten is probably advisable. Values of 3 or 4 work well for initial trouble shooting and testing. 
#' @param report_freq Divisor specifying report frequency during optimization; will print updates every \code{report_freq} accepted calls, default 10.
#' @param maxeval Maximum number of objective function evaluations for L-BFGS, default 500.
#' @param maxHam Maximum Hamming distance for misreads to be corrected (max correctable Hamming distance), default NULL sets it to one less than the minimum Hamming distance between codebook entries.
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
    weights              = list(
      blank      = 1.0, 
      prior      = 0.0, 
      dispersion = 0.0, 
      erc        = 0.0
    )
  ) {
    cat("\nBenchmarking misread QC with Bernoulli simulations, using L-BFGS (nlopt)")
   
    # Fill in any weight the caller's list omits with its default (same rationale as misread.qc())
    weights <- modifyList(
      list(blank = 1.0, prior = 0.0, dispersion = 0.0, erc = 0.0),
      weights
    )
   
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
      as.double(weights$blank),
      as.double(weights$prior),
      as.double(weights$dispersion),
      as.double(weights$erc),
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
