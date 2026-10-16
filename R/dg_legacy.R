
# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2026

# #########################################################################################################
# Verbatim port of the original "dichotomized Gaussian" (DG) spot-simulation method -- the simulation-based
# approach used to produce this project's earlier results -- from the external repo/folder
# paper_calling_rate_sim/MERFISH_calling_rate_sim/main.R (paired with src/dg_legacy.cpp, the ported
# fast_functions.cpp). The exported entry point is misread.qc.DG.legacy() (renamed from the original
# misread.qc.DG() when the C++-backed, actively-developed misread.qc.DG()/mQC_DG() was added to
# R/main.R / src/main.cpp). This file exists only so the original results can in principle be reproduced;
# it is not the package's primary (recommended) method. Nothing in main.R's numerical logic is shared with
# or calls into this file -- the one exception is that misread.qc.DG()'s thin R wrapper reuses this file's
# private .DG_plot_* ggplot helpers for its own post-hoc plotting, since both files are part of the same
# package namespace; those helpers themselves are untouched.
#
# Every internal helper below is a straight rename of its original counterpart to a ".DG_"-prefixed name
# (so it stays unexported under the package's `exportPattern("^[[:alpha:]]+")` and can never collide with
# anything else in flipFISH); every function BODY is otherwise unchanged. The only two substantive
# deviations from verbatim, both confined to the single exported entry point below, are:
#   1. Input data: the original top-level driver (`run_MCMCSA`) read summary-stats CSVs from disk by
#      mouse_id/path. misread.qc.DG.legacy() instead takes STdata/codebook in-memory, in the same shape
#      misread.qc() already uses, and splits it locally -- no change to any downstream numerical logic.
#   2. Parallelism: the original relied on a script-level global (`core_num`) read inside `runparallel()`.
#      Here `n_forks` is threaded through explicitly as a function argument instead.
# Note that maxHam has also been exposed as an argument to misread.qc.DG.legacy with a default value of 4
# matching the original hard-coded value.
# #########################################################################################################

# ----- style constants (originally script-level globals used by the plotting helpers) --------------------

.DG_title_size           <- 18
.DG_axis_size            <- 14
.DG_legend_size          <- 14
.DG_dot_color             <- "black"
.DG_filter_color          <- "magenta"
.DG_filter_color_obs      <- "magenta4"
.DG_filter_color2         <- "cyan"
.DG_filter_color2_obs     <- "cyan4"
.DG_filter_style          <- "solid"
.DG_retained_color        <- "skyblue1"
.DG_retained_color_obs    <- "skyblue4"
.DG_filtered_color        <- "firebrick1"
.DG_filtered_color_obs    <- "firebrick4"
.DG_blanks_color          <- "gray"
.DG_blanks_color_obs      <- "gray20"

# Helper function to check forks (mirrors check_forks() in main.R; duplicated rather than shared so this
# legacy file has no dependency on main.R)
.DG_check_forks <- function(
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

# For lazy parallelization (originally used a script-level global `core_num` via lexical scope; here
# n_forks is passed explicitly instead -- the mclapply() call itself is unchanged)
.DG_runparallel <- function(input, funct, n_forks = 1, ...) {
    results <- parallel::mclapply(
      X = input, FUN = funct, ...,
      mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = n_forks,
      mc.cleanup = TRUE, mc.allow.recursive = TRUE )
    return(results)
  }

# 3. Functions to generate simulated data ##############################################################

# Helper function to set true counts for each barcode
.DG_set_true_spot_info <- function(
    gene_kernel_rates,
    gamma_variance,
    blank_mask,
    N_cells,
    N_barcodes,
    N_blanks,
    N_genes
  ) {
   
    # Stipulate true transcript counts per cell, using over-dispersed Poisson distribution
    # ... translate these values into the standard gamma-distribution parameters (rate and shape)
    gamma_rt <- gene_kernel_rates/gamma_variance
    gamma_sp <- gene_kernel_rates^2/gamma_variance
    # ... set transcription rates from the gamma distribution using the kernel rate and gamma variance
    transcription_rates <- matrix(NA, nrow = N_cells, ncol = N_genes)
    for (ci in 1:N_cells) {
      transcription_rates[ci,] <- rgamma(N_genes, shape = gamma_sp, rate = gamma_rt)
    }
    # ... convert this transcription rate into an integer count by drawing from a Poisson distribution of that expected value
    barcode_counts_true <- matrix(NA, nrow = N_cells, ncol = N_barcodes)
    for (ci in 1:N_cells) {
      barcode_counts_true[ci,!blank_mask] <- rpois(N_genes, lambda = transcription_rates[ci,])
    }
   
    # Set all true blank counts to zero
    blank_idx <- which(blank_mask)
    for (bci in 1:N_blanks) {
      barcode_counts_true[,blank_idx[bci]] <- 0
    }
   
    # Collect true barcode counts for each cell
    barcodes_by_cell_true <- list()
    length(barcodes_by_cell_true) <- N_cells
    for (ci in 1:N_cells) {
      barcodes_by_cell_true[[ci]] <- rep(1:ncol(barcode_counts_true), times = barcode_counts_true[ci,])
      # ... this is n_spots, or sum(barcode_counts_true[ci,])
    }
   
    # Return results:
    #   barcode_counts_true: N_cells x N_barcodes matrix giving true counts of each barcode in each cell
    #   barcodes_by_cell_true: list of length N_cells, each element is a vector of length n_spots (varies by cell), giving the barcode index for each spot in that cell
    out <- list(
      barcode_counts_true   = barcode_counts_true,
      barcodes_by_cell_true = barcodes_by_cell_true
    )
   
  }

# Function for making codebook
.DG_make_codebook_generation_matrix <- function(
    codebook
  ) {
    # ... for min and max in runif, the exact range doesn't matter, as flip rates are controlled via the variance
    codebook_generation_matrix <- codebook
    for (i in 1:nrow(codebook_generation_matrix)) {
      for (j in 1:ncol(codebook_generation_matrix)) {
        if (codebook[i,j] == 0) {
          codebook_generation_matrix[i,j] <- runif(1, min = -1, max = 0)
        } else {
          codebook_generation_matrix[i,j] <- runif(1, min = 1e-12, max = 1)
        }
      }
    }
    codebook_generation_matrix <- as.matrix(codebook_generation_matrix)
    return(codebook_generation_matrix)
  }

# Function to set the correlation matrix for luminance noise
.DG_set_luminance_noise_correlation_matrix <- function(
    noise_correlation_bias = 0,
    noise_correlation_sd   = 1,
    lum_noise_correlations = NULL,
    N_bits,
    N_corrs
  ) {
    # ... randomly select correlations, if needed
    if (is.null(lum_noise_correlations)) {
      lum_noise_correlations <- rnorm(N_corrs, mean = noise_correlation_bias, sd = noise_correlation_sd)
    }
    # ... clamp these (nearPD would force this anyway)
    lum_noise_correlations[lum_noise_correlations < -1] <- -1
    lum_noise_correlations[lum_noise_correlations > 1]  <- 1
    # ... distribute into correlation matrix
    Sigma                   <- matrix(NA, nrow = N_bits, ncol = N_bits)
    Sigma[upper.tri(Sigma)] <- lum_noise_correlations
    # ... make diagonal 1
    diag(Sigma) <- 1
    # ... fill bottom triangle via symmetry
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
    # ... enforce positive definiteness
    Sigma <- Matrix::nearPD(Sigma, corr = TRUE, keepDiag = TRUE)$mat
    return(Sigma)
  }

# Function to set noise by flip rates
.DG_set_noise_by_flip_rates <- function(
    codebook_generation_matrix,
    flip_rate_10,
    flip_rate_01
  ) {
   
    # Initialize matrix to hold noise values (sd for each bit, for each barcode)
    codebook_dim   <- dim(codebook_generation_matrix)
    N_bits         <- codebook_dim[2]
    N_barcodes     <- codebook_dim[1]
    codebook_noise <- matrix(NA, nrow = N_barcodes, ncol = N_bits)
   
    # Check that flip rates are in the right form
    length_flip_rate_10 <- length(flip_rate_10)
    if (length_flip_rate_10 != N_bits) {
      if (length_flip_rate_10 == 1) {
        flip_rate_10 <- rep(flip_rate_10, N_bits)
      } else {
        stop("flip_rate_10 must be a single value or a vector of length equal to the number of bits in the codebook.")
      }
    }
    length_flip_rate_01 <- length(flip_rate_01)
    if (length_flip_rate_01 != N_bits) {
      if (length_flip_rate_01 == 1) {
        flip_rate_01 <- rep(flip_rate_01, N_bits)
      } else {
        stop("flip_rate_01 must be a single value or a vector of length equal to the number of bits in the codebook.")
      }
    }
   
    for (bi in 1:N_bits) {
      # For bit bi, grab the seeds (expected luminance value) for each barcode
      these_seeds                 <- codebook_generation_matrix[,bi]
      # Make mask for the positive bits (1)
      pos_mask                    <- these_seeds > 0
      # Make mask for the negative bits (0)
      neg_mask                    <- these_seeds < 0
      # Find luminance value such that flip_rate_10 of a standard normal distribution would fall below this value
      inv_flip_rate_10            <- qnorm(flip_rate_10[bi])
      # "denormalize" to find the standard deviation of the distribution with this same qnorm, but the seeded mean
      codebook_noise[pos_mask,bi] <- -these_seeds[pos_mask] / inv_flip_rate_10
      # Same as above, but for flip_rate_01
      inv_flip_rate_01            <- qnorm(1 - flip_rate_01[bi])
      codebook_noise[neg_mask,bi] <- -these_seeds[neg_mask] / inv_flip_rate_01
    }
   
    return(codebook_noise)
   
  }

# Function to extract luminance noise from spots
.DG_find_noise <- function(
    ci,
    spots,
    bc_counts,
    codebook_gen_mat,
    N_bits
  ) {
    spots_ci     <- spots[[ci]]
    N_barcodes   <- ncol(bc_counts)
    total_rows   <- sum(bc_counts[ci,])
    noise_ci     <- matrix(NA, nrow = total_rows, ncol = N_bits)
    # Precompute start and end indices per barcode
    bc_counts_ci <- bc_counts[ci, ]
    end_rows     <- cumsum(bc_counts_ci)
    start_rows   <- end_rows - bc_counts_ci + 1
    for (bci in seq_len(N_barcodes)) {
      s <- start_rows[bci]
      e <- end_rows[bci]
      if (e >= s) {
        noise_ci[s:e, ] <- spots_ci[s:e, , drop = FALSE] - matrix(codebook_gen_mat[bci, ], nrow = e - s + 1, ncol = N_bits, byrow = TRUE)
      }
    }
    return(noise_ci)
  }

# Function to generate simulated spots
.DG_generate_simulated_spots <- function(
    ci,
    noise,
    bc_counts,
    codebook_gen_mat,
    corr_mat,
    N_bits,
    decode_and_label,
    codebook  = NULL,
    maxHam    = NULL,
    threshold = 0
  ) {
    N_barcodes <- ncol(bc_counts)
   
    use_Rcpp <- TRUE
    if (use_Rcpp) {
      spots_ci <- generate_spots_ci(
        as.matrix(bc_counts) * 1L,
        as.matrix(codebook_gen_mat),
        as.matrix(noise),
        as.matrix(corr_mat),
        ci - 1,
        threshold,
        decode_and_label
      )
    } else {
      # ... initialize matrix for spots
      spots_ci <- matrix(NA, nrow = sum(bc_counts[ci,]), ncol = N_bits)
      for (bci in 1:N_barcodes) {
        # ... for each barcode, use true counts to find its rows in spots_ci
        end_row   <- sum(bc_counts[ci,1:bci])
        start_row <- end_row - bc_counts[ci,bci] + 1
        # ... if there are rows to fill, generate spots
        if (end_row >= start_row) {
          # ... grab the generation kernel values for this barcode
          this_bc      <- codebook_gen_mat[bci,]
          # ... set the sd for the noise for each bit and put into matrix form
          noise_matrix <- diag(noise[bci,])
          # ... turn correlation matrix into covariance matrix with proper sd for the noise causing the stipulated flip rates
          Sigma_noised <- noise_matrix %*% corr_mat %*% noise_matrix
          # ... generate noised kernel
          spot_seeds   <- MASS::mvrnorm(n = bc_counts[ci,bci], mu = this_bc, Sigma = Sigma_noised)
          # ... save
          if (decode_and_label) spots_ci[start_row:end_row,] <- matrix(as.numeric(spot_seeds > threshold), nrow = bc_counts[ci,bci])
          else spots_ci[start_row:end_row,]                  <- matrix(spot_seeds, nrow = bc_counts[ci,bci])
        }
      }
    }
    if (decode_and_label) {
      # Call C++ function which finds corrections
      spot_labels <- find_spot_labels_fast(
        as.matrix(codebook) * 1L,
        as.matrix(spots_ci) * 1L,
        maxHam * 1L
      )
      # Returns a list of two integer vectors, the first corrected labels, the second uncorrected labels
      return(spot_labels)
    } else {
      return(spots_ci)
    }
   
  }

# Function to generate RNA molecules
.DG_generate_RNA <- function(
    ci,
    bc_counts,
    codebook,
    N_bits
  ) {
    # The function .DG_generate_simulated_spots generates the observed luminance sequences, while
    #   this function generates the corresponding "ground truth" RNA molecules generating those sequences
   
    # Get the counts for the current cell
    counts          <- bc_counts[ci, ]
    total_RNA       <- sum(counts)
    # Repeat the barcode indices according to their counts
    barcode_indices <- rep(seq_len(nrow(codebook)), times = counts)
    # Construct the RNA matrix by indexing the codebook
    RNA_ci          <- matrix(as.numeric(t(codebook[barcode_indices, ])), nrow = total_RNA, ncol = N_bits, byrow = TRUE)
   
    return(RNA_ci)
    
  }

# 4. Functions to decode and correct observed spots #####################################################

# Function to decode spots
.DG_decode_spots <- function(
    ci,
    spots,
    threshold = 0
  ) {
    spots_ci <- spots[[ci]]
    spots_ci <- matrix(as.numeric(spots_ci > threshold), nrow = nrow(spots_ci), ncol = ncol(spots_ci))
    return(spots_ci)
  }

# Function to find percentage of spots decoded
.DG_find_decoded_spot_rate <- function(
    spot_labels
  ) {
    N_cells       <- length(spot_labels)
    total_spots   <- 0
    decoded_spots <- 0
    for (ci in 1:N_cells) {
      total_spots   <- total_spots + length(spot_labels[[ci]])
      decoded_spots <- decoded_spots + sum(!is.na(spot_labels[[ci]]))
    }
    return(decoded_spots / total_spots)
  }

# Extract barcode counts from the labels
.DG_find_barcode_counts <- function(
    spot_bc_labels_list,
    N_barcodes
  ) {
    barcode_counts_est <- t(sapply(
      spot_bc_labels_list,
      function(labels) tabulate(labels, nbins = N_barcodes)
    ))
    return(barcode_counts_est)
  }

# Function to correct generated spots
.DG_apply_spot_correction <- function(
    ci,
    spots,
    codebook,
    maxHam
  ) {
    # Get matrix of spots
    spots_ci <- spots[[ci]]
    # Call C++ function which finds corrections
    spot_labels <- find_spot_labels_fast(
      as.matrix(codebook) * 1L,
      as.matrix(spots_ci) * 1L,
      maxHam * 1L
      )
    # Returns a list of two integer vectors, the first corrected labels, the second uncorrected labels
    return(spot_labels)
  }

# Gathers the unique distances corrected for a given cell
.DG_check_correctable_distances <- function(
    spots,
    codebook,
    maxHam
  ) {
   
    N_cells <- length(spots)
    distances_found <- rep(FALSE, maxHam + 1)
   
    for (ci in 1:N_cells) {
     
      spots_ci <- spots[[ci]]
      n_spots <- nrow(spots_ci)
      for (si in 1:n_spots) {
        spot <- spots_ci[si,]
        # ... find distance from this spot to all barcodes in the codebook
        hamming_distances <- codebook != matrix(spot, nrow = nrow(codebook), ncol = ncol(codebook), byrow = TRUE)
        hamming_distances <- rowSums(hamming_distances)
        # ... apply correction (if possible) to the observed spot
        min_HD <- min(hamming_distances)
        if (min_HD <= maxHam) {
          min_HD_mask <- hamming_distances == min_HD
          if (sum(min_HD_mask) == 1) {
            distances_found[min_HD + 1] <- TRUE
          }
        }
        if (all(distances_found)) break
      }
     
      if (all(distances_found)) break
     
    }
   
    distances_found <- c(0:maxHam)[distances_found]
    return(distances_found)
   
  }

# Calculate flip rates per cell (without correction)
.DG_check_flip_rates <- function(
    RNA,
    spots,
    verbose = TRUE
  ) {
   
    # Find flip rates and take mean
    N_cells   <- length(RNA)
    pos_flips <- 0
    pos_bits  <- 0
    neg_flips <- 0
    neg_bits  <- 0
    for (ci in 1:N_cells) {
      RNA_ci   <- c(RNA[[ci]])
      spots_ci <- c(spots[[ci]])
      pos_mask <- RNA_ci > 0
      neg_mask <- RNA_ci == 0
      # ... number of 1s that flipped to 0 and total number of 1s
      pos_flips <- pos_flips + sum(spots_ci[pos_mask] != 1)
      pos_bits  <- pos_bits + sum(pos_mask)
      # ... number of 0s that flipped to 1 and total number of 0s
      neg_flips <- neg_flips + sum(spots_ci[neg_mask] != 0)
      neg_bits  <- neg_bits + sum(neg_mask)
    }
    flip_rate_10_mean <- pos_flips / pos_bits
    flip_rate_01_mean <- neg_flips / neg_bits
    if (verbose) {
      cat("\n\tFlip rate 1 -> 0:", flip_rate_10_mean)
      cat("\n\tFlip rate 0 -> 1:", flip_rate_01_mean)
    }
   
    out <- c(flip_rate_10_mean, flip_rate_01_mean)
    names(out) <- c("flip_rate_10", "flip_rate_01")
    return(out)
   
  }

# 5. Functions to run and analyze simulation #############################################################

# Function to run individual simulation
.DG_run_simulation <- function(
    codebook,
    summary_stats_genes,
    summary_stats_blanks,
    flip_rate_10,                           # Probability of a 1 being flipped to a 0, either scalar or vector of length equal to number of bits
    flip_rate_01,                           # Probability of a 0 being flipped to a 1, either scalar or vector of length equal to number of bits
    maxHam,                                 # For the automatic error correction, assumes an "extended" Hamming code with min 4 bits between all barcodes
    bit_noise_correlations_bias = 0,        # For correlating noise in the bit reads
    bit_noise_correlations_sd   = 1,        # For correlating noise in the bit reads
    bit_lum_noise_correlations  = NULL,     # For stipulating correlations between bits (bit bias and sd will be ignored)
    N_cells                     = 1e3,      # Number of cells to simulate
    blank_weight                = 1,
    maintain_gene_identity      = TRUE,
    return_flip_rates           = FALSE,
    n_forks                     = 1,
    verbose                     = FALSE
  ) {
   
    # Extract summary stats
    observed_gene_rates     <- summary_stats_genes$rates
    observed_gene_variance  <- summary_stats_genes$variance
    observed_blank_rates    <- summary_stats_blanks$rates
    observed_blank_variance <- summary_stats_blanks$variance
   
    # Extract simulation parameters
    N_bits              <- ncol(codebook)
    N_barcodes          <- nrow(codebook)
    N_blanks            <- length(observed_blank_rates)
    N_genes             <- length(observed_gene_rates)
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
      cat("\n\tMax correctable Hamming distance:", maxHam)
    }
   
    # Make generation matrix for codebook, representing the expected luminance values for each bit in each barcode
    # ... regenerated with random (but properly valenced) numbers each simulation
    if (verbose) cat("\nMaking codebook generation matrix")
    codebook_generation_matrix <- .DG_make_codebook_generation_matrix(codebook)
   
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
    Sigma <- .DG_set_luminance_noise_correlation_matrix(
      noise_correlation_bias = bit_noise_correlations_bias,
      noise_correlation_sd   = bit_noise_correlations_sd,
      lum_noise_correlations = bit_lum_noise_correlations,
      N_bits                 = N_bits,
      N_corrs                = N_corrs
    )
    
    # Predict decoding rate
    expected_flip_rate     <- (mean(flip_rate_10)*mean_Hamming_weight + mean(flip_rate_01)*(N_bits - mean_Hamming_weight)) / N_bits
    expected_decoding_rate <- ppois(q = maxHam/2, lambda = expected_flip_rate * N_bits)
    if (verbose) cat("\nPredicted decoding rate:", expected_decoding_rate)
   
    # Set true counts
    # ... set randomly each simulation, based around the summary statistics of the observed data
    if (verbose) cat("\nSetting true counts for each gene and blank")
    true_spot_info <- .DG_set_true_spot_info(
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
    codebook_noise <- .DG_set_noise_by_flip_rates(
      codebook_generation_matrix = codebook_generation_matrix,
      flip_rate_10               = flip_rate_10,
      flip_rate_01               = flip_rate_01
    )
   
    # Generate simulated spots
    if (verbose) cat("\nGenerating simulated spots")
    spots <- .DG_runparallel(
      1:N_cells,
      .DG_generate_simulated_spots,
      n_forks          = n_forks,
      noise            = codebook_noise,
      bc_counts        = true_spot_info$barcode_counts_true,
      codebook_gen_mat = codebook_generation_matrix,
      corr_mat         = Sigma,
      N_bits           = N_bits,
      decode_and_label = !verbose && !return_flip_rates,
      codebook         = codebook,
      maxHam           = maxHam
    )
   
    # ... initialize flip rates
    # ... computationally, very expensive to check, only check if verbose
    flip_rates <- c(NA, NA)
    names(flip_rates) <- c("flip_rate_10", "flip_rate_01")
    
    if (verbose) {
     
      # Check luminance noise correlations
      cat("\nChecking luminance noise correlations by bit")
      # ... find noise by bit by cell
      luminance_noise <- .DG_runparallel(
        1:N_cells,
        .DG_find_noise,
        n_forks          = n_forks,
        spots            = spots,
        bc_counts        = true_spot_info$barcode_counts_true,
        codebook_gen_mat = codebook_generation_matrix,
        N_bits           = N_bits
      )
      # ... collapse into single matrix of noise values (rows as spots, columns as bits)
      luminance_noise <- do.call(rbind, luminance_noise)
      # ... calculate correlations in noise between bits
      corrs_all       <- matrix(NA, nrow = N_corrs, ncol = 2)
      idx             <- 0
      for (bi in 1:(N_bits-1)) {
        for (bj in (bi+1):N_bits) {
          idx                  <- idx + 1
          noise_correlation <- cor(
            luminance_noise[,bi],
            luminance_noise[,bj]
          )
          expected_correlation <- Sigma[bi,bj]
          corrs_all[idx,]      <- c(noise_correlation, expected_correlation)
        }
      }
      # ... print results
      cat("\nCorrelation between observed and expected luminance-noise correlations:",
          cor(corrs_all[,1], corrs_all[,2], use = "pairwise.complete.obs"))
     
      # Decode spots
      cat("\nDecoding spots")
      spots <- .DG_runparallel(
        1:N_cells,
        .DG_decode_spots,
        n_forks   = n_forks,
        spots     = spots,
        threshold = 0
      )
     
      # Check unique correctable error sizes
      cat("\nChecking unique correctable error sizes")
      corrected_distances <- .DG_check_correctable_distances(
        spots    = spots,
        codebook = codebook,
        maxHam   = maxHam
      )
      cat("\nUnique correctable error sizes:",
          paste(corrected_distances, collapse = ", "))
     
      # Find the ground truth behind the simulated spots
      cat("\nGenerating simulated mRNA molecules responsible for those spots (ground truth)")
      RNA <- .DG_runparallel(
        1:N_cells,
        .DG_generate_RNA,
        n_forks   = n_forks,
        bc_counts = true_spot_info$barcode_counts_true,
        codebook  = codebook,
        N_bits    = N_bits
      )
     
      # Check flip rates
      cat("\nChecking bit-flip rates")
      flip_rates <- .DG_check_flip_rates(
        RNA,
        spots,
        verbose = verbose
      )
     
      # Find uncorrected bc and apply read corrections
      if (verbose) cat("\nFinding uncorrected barcode labels and applying read corrections")
      spots <- .DG_runparallel(
        1:N_cells,
        .DG_apply_spot_correction,
        n_forks  = n_forks,
        spots    = spots,
        codebook = codebook,
        maxHam   = maxHam
      )
     
    } else if (return_flip_rates) {
     
      # Decode spots
      spots <- .DG_runparallel(
        1:N_cells,
        .DG_decode_spots,
        n_forks   = n_forks,
        spots     = spots,
        threshold = 0
      )
     
      # Find the ground truth behind the simulated spots
      RNA <- .DG_runparallel(
        1:N_cells,
        .DG_generate_RNA,
        n_forks   = n_forks,
        bc_counts = true_spot_info$barcode_counts_true,
        codebook  = codebook,
        N_bits    = N_bits
      )
     
      # Check flip rates
      flip_rates <- .DG_check_flip_rates(
        RNA,
        spots,
        verbose = verbose
      )
     
      # Find uncorrected bc and apply read corrections
      spots <- .DG_runparallel(
        1:N_cells,
        .DG_apply_spot_correction,
        n_forks  = n_forks,
        spots    = spots,
        codebook = codebook,
        maxHam   = maxHam
      )
      
    }
    
    spots_corrected_bc_labels    <- lapply(spots, `[[`, "spot_labels")
    spots_bc_labels              <- lapply(spots, `[[`, "spot_labels_uncorrected")
    spot_decoding_rate           <- .DG_find_decoded_spot_rate(spots_bc_labels)
    spot_decoding_rate_corrected <- .DG_find_decoded_spot_rate(spots_corrected_bc_labels)
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
    barcode_counts_est                        <- .DG_find_barcode_counts(spots_bc_labels, N_barcodes)
    barcode_counts_est_corrected              <- .DG_find_barcode_counts(spots_corrected_bc_labels, N_barcodes)
    # ... extract results
    total_barcode_counts_est                  <- colSums(barcode_counts_est, na.rm = TRUE)
    total_barcode_counts_est_genes            <- total_barcode_counts_est[!blank_mask]
    total_barcode_counts_est_corrected        <- colSums(barcode_counts_est_corrected, na.rm = TRUE)
    total_barcode_counts_est_corrected_genes  <- total_barcode_counts_est_corrected[!blank_mask]
    total_barcode_counts_est_corrected_blanks <- total_barcode_counts_est_corrected[blank_mask]
    max_blank_count                           <- max(total_barcode_counts_est_corrected_blanks)
    sim_gene_order                            <- order(total_barcode_counts_est_corrected_genes, decreasing = TRUE)
    sim_blank_order                           <- order(total_barcode_counts_est_corrected_blanks, decreasing = TRUE)
    
    # ... while we're at it, make sorted vector for the real data
    obs_counts_genes                          <- observed_gene_rates * N_cells
    obs_counts_blanks                         <- observed_blank_rates * N_cells
    obs_counts                                <- rep(NA, N_barcodes)
    obs_counts[!blank_mask]                   <- obs_counts_genes
    obs_counts[blank_mask]                    <- obs_counts_blanks
    max_blank_count_obs                       <- max(obs_counts_blanks)
    obs_gene_order                            <- order(obs_counts_genes, decreasing = TRUE)
    obs_blank_order                           <- order(obs_counts_blanks, decreasing = TRUE)
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
    weights     <- c(rep(gene_weight, N_genes), rep(blank_weight, N_blanks))
    sim_se      <- (log(counts_sorted + 1) - log(obs_counts_sorted + 1))^2
    sim_mse     <- mean(sim_se * weights)
    if (verbose) cat("\nWeighted mean squared log-error:", sim_mse)
    
    # Compute confidence ratio for barcode counts
    if (verbose) cat("\nComputing confidence ratio")
    CR <- total_barcode_counts_est / total_barcode_counts_est_corrected
    CR[total_barcode_counts_est_corrected == 0] <- 0
    if (any(is.na(CR) | is.infinite(CR) | is.nan(CR) | isTRUE(CR <= 0))) {
      stop("Confidence ratio contains invalid values (NA, Inf, NaN, or <= 0).")
    }
    # ... extract results
    CR_genes     <- CR[!blank_mask]
    CR_blanks    <- CR[blank_mask]
    max_blank_CR <- max(CR_blanks)
    CR_sorted    <- c(
      CR_genes[order(CR_genes, decreasing = TRUE)],
      CR_blanks[order(CR_blanks, decreasing = TRUE)]
    )
   
    # Find relationship between CR and PPV
    # ... fit linear model with intercept = 0
    CR_PPV_lm       <- lm(PPV_corrected_decoding ~ CR_genes + 0)
    # ... extract slope
    CR_PPV_lm_slope <- coef(CR_PPV_lm)[["CR_genes"]]
    # ... compute predictions and MSE
    predictions     <- CR_PPV_lm_slope * CR_genes
    mse             <- mean((PPV_corrected_decoding - predictions)^2)
    if (verbose) {
      cat("\nLinear model of PPV ~ CR + 0:")
      cat("\n\tSlope:", CR_PPV_lm_slope)
      cat("\n\tMSE:", mse)
    }
   
    # Extract cutoff stats
    # ... assume that PPV is normally distributed for genes below the max-count blank cutoff
    retained_gene_mask      <- total_barcode_counts_est_corrected_genes > max_blank_count
    PPV_mean_filtered       <- mean(PPV_corrected_decoding[!retained_gene_mask])
    PPV_sd_filtered         <- sd(PPV_corrected_decoding[!retained_gene_mask])
    # ... assume that PPV is exponentially distributed for genes above the max-count blank cutoff
    PPV_mean_retained       <- mean(PPV_corrected_decoding[retained_gene_mask])
    # ... get percent of genes filtered by CR cutoff
    ratio_genes_CR_filtered <- sum(CR_genes <= max_blank_CR) / length(CR_genes)
   
    # Gather results
    bit_noise_correlations_bias_out <- bit_noise_correlations_bias
    bit_noise_correlations_sd_out   <- bit_noise_correlations_sd
    if (!is.null(bit_lum_noise_correlations)) {
      bit_noise_correlations_bias_out <- mean(bit_lum_noise_correlations)
      bit_noise_correlations_sd_out   <- sd(bit_lum_noise_correlations)
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
      Count           = total_barcode_counts_est_genes,
      Count_corrected = total_barcode_counts_est_corrected_genes,
      CR              = CR_genes,
      PPV             = PPV_corrected_decoding
    )
    rownames(gene_summary) <- rownames(codebook)[!blank_mask]
   
    results <- list(
      gene_summary  = gene_summary,
      counts_sorted = data.frame(sim = counts_sorted, obs = obs_counts_sorted),
      sim_summary   = sim_summary
    )
   
    if (verbose) {
     
      blank_line_width      <- 1.0
      observed_gene_counts  <- summary_stats_genes$counts
      observed_blank_counts <- summary_stats_blanks$counts
     
      # barcode count distributions
      df_obs_est  <- data.frame(x = log(c(observed_gene_counts, observed_blank_counts)+1))
      df_sim_est  <- data.frame(x = log(colSums(barcode_counts_est_corrected)+1))
      df_sim_true <- data.frame(x = log(colSums(true_spot_info$barcode_counts_true)+1))
     
      cat("\nMean spots per cell:",
          "\n\tObserved, decoded:",  mean(c(observed_gene_rates, observed_blank_rates)),
          "\n\tSimulated, decoded:", mean(colSums(barcode_counts_est_corrected)/N_cells),
          "\n\tSimulated, true:",    mean(colSums(true_spot_info$barcode_counts_true)/N_cells))
     
      plot_obs_est <- ggplot(df_obs_est, aes(x)) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black") +
        labs(title = "Real decoded counts", x = "Barcode count", y = NULL) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size*0.8),
          axis.title      = element_text(size = .DG_axis_size),
          axis.text       = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      plot_sim_est <- ggplot(df_sim_est, aes(x)) +
        geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
        labs(title = "Simulated decoded counts", x = "Barcode count", y = NULL) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size*0.8),
          axis.title      = element_text(size = .DG_axis_size),
          axis.text       = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      plot_sim_true <- ggplot(df_sim_true, aes(x)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        labs(title = "Simulated ground-truth counts", x = "Barcode count", y = NULL) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size*0.8),
          axis.title      = element_text(size = .DG_axis_size),
          axis.text       = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      gridExtra::grid.arrange(
        plot_obs_est, plot_sim_est, plot_sim_true,
        nrow = 1,
        top = grid::textGrob(
          "Barcode Count Frequency",
          gp = grid::gpar(fontsize = .DG_title_size)
        )
      )
     
      bc_type_colors <- c(
        "Retained" = .DG_retained_color, "Filtered" = .DG_filtered_color, "Blanks" = .DG_blanks_color,
        "Retained_obs" = .DG_retained_color_obs, "Filtered_obs" = .DG_filtered_color_obs, "Blanks_obs" = .DG_blanks_color_obs
      )
     
      # CR filter
      bc_type                                          <- rep("Retained", length(CR_sorted))
      bc_type[CR_sorted <= max_blank_CR]               <- "Filtered"
      bc_type[seq_along(CR_sorted) > length(CR_genes)] <- "Blanks"
      df_CR_sorted <- data.frame(
        index   = seq_along(CR_sorted),
        CR      = CR_sorted,
        bc_type = bc_type
      )
      df_CR_sorted$CR[df_CR_sorted$CR == 0] <- min(df_CR_sorted$CR[df_CR_sorted$CR > 0]) * 0.01
      df_CR_sorted$bc_type                  <- factor(df_CR_sorted$bc_type, levels = c("Retained", "Filtered", "Blanks"))
      plot_CR_sorted <- ggplot(df_CR_sorted, aes(x = index, y = CR, color = bc_type)) +
        geom_point(size = 3) +
        scale_y_log10() +  # Log scale for better visibility
        scale_color_manual(values = bc_type_colors) +
        geom_vline(
          xintercept = mean(c(length(CR_genes), length(CR_genes)+1)),
          color      = "black",
          linetype   = "dashed",
          linewidth  = 0.5
        ) +
        geom_vline(
          xintercept = mean(c(min(which(CR_sorted <= max_blank_CR)), min(which(CR_sorted <= max_blank_CR)) + 1)),
          color      = .DG_filter_color,
          linetype   = .DG_filter_style,
          linewidth  = blank_line_width
        ) +
        labs(title = "Max-CR filter", x = "Barcodes sorted by CR", y = "Confidence ratio (CR)", color = "Barcode status") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size),
          axis.title      = element_text(size = .DG_axis_size),
          axis.text.y     = element_text(size = .DG_axis_size),
          axis.text.x     = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      # Count filter
      bc_type                                                                              <- rep("Retained", length(counts_sorted))
      bc_type[counts_sorted <= max_blank_count]                                            <- "Filtered"
      bc_type[seq_along(counts_sorted) > length(total_barcode_counts_est_corrected_genes)] <- "Blanks"
      bc_type_obs                                                                          <- rep("Retained_obs", length(obs_counts_sorted))
      bc_type_obs[obs_counts_sorted <= max_blank_count_obs]                                <- "Filtered_obs"
      bc_type_obs[seq_along(obs_counts_sorted) > length(obs_counts_genes)]                 <- "Blanks_obs"
      df_counts_sorted <- data.frame(
        index     = c(seq_along(counts_sorted), seq_along(obs_counts_sorted)),
        Count     = c(counts_sorted, obs_counts_sorted),
        bc_type   = c(bc_type, bc_type_obs),
        data_type = c(rep("Simulated", length(counts_sorted)), rep("Observed", length(obs_counts_sorted))),
        size      = c(rep(3, length(counts_sorted)), rep(1.5, length(obs_counts_sorted)))
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
          color      = .DG_filter_color_obs,
          linetype   = .DG_filter_style,
          linewidth  = blank_line_width/2
        ) +
        geom_vline(
          xintercept = mean(c(length(total_barcode_counts_est_corrected_genes), length(total_barcode_counts_est_corrected_genes)+1)),
          color      = "black",
          linetype   = "dashed",
          linewidth  = 0.5
        ) +
        geom_vline(
          xintercept = mean(c(min(which(counts_sorted <= max_blank_count)), min(which(counts_sorted <= max_blank_count)) + 1)),
          color      = .DG_filter_color,
          linetype   = .DG_filter_style,
          linewidth  = blank_line_width
        ) +
        labs(title = "Max-count filter", x = "Barcodes sorted by count", y = "Total counts per barcode", color = "Barcode status") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size),
          axis.title      = element_text(size = .DG_axis_size),
          axis.text.y     = element_text(size = .DG_axis_size),
          axis.text.x     = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      print(plot_CR_sorted)
      print(plot_counts_sorted)
     
      # CR vs PPV (genes only)
      ret                          <- rep("no", length(CR_genes))
      ret[CR_genes > max_blank_CR] <- "yes"
      df_CR_filter                 <- data.frame(
        CR       = CR_genes,
        PPV      = PPV_corrected_decoding,
        Retained = ret
      )
      plot_CR_filter <- ggplot(df_CR_filter, aes(x = CR, y = PPV, color = Retained)) +
        geom_point() +
        scale_color_manual(values = c("no" = "gray", "yes" = .DG_dot_color)) +
        geom_vline(
          xintercept = max_blank_CR,
          color      = .DG_filter_color,
          linetype   = .DG_filter_style,
          linewidth  = blank_line_width
        ) +
        ylim(c(0, 1)) +  # Limit y-axis to [0, 1]
        xlim(c(0, 1)) +
        labs(title = "Max blank CR", x = "Confidence ratio (CR)", y = NULL) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size*0.8),
          axis.title.     = element_text(size = .DG_axis_size),
          axis.text       = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      # Total Barcode Counts vs PPV (genes only)
      ret                                                             <- rep("no", length(total_barcode_counts_est_corrected_genes))
      ret[total_barcode_counts_est_corrected_genes > max_blank_count] <- "yes"
      df_count_filter <- data.frame(
        Count    = total_barcode_counts_est_corrected_genes,
        PPV      = PPV_corrected_decoding,
        Retained = ret
      )
      plot_count_filter <- ggplot(df_count_filter, aes(x = Count, y = PPV, color = Retained)) +
        geom_point() +
        scale_x_log10() +  # Log scale for better visibility
        scale_color_manual(values = c("no" = "gray", "yes" = .DG_dot_color)) +
        geom_vline(
          xintercept = max_blank_count,
          color      =  .DG_filter_color,
          linetype   = .DG_filter_style,
          linewidth  = blank_line_width
        ) +
        labs(title = "Max blank count", x = "Total counts per barcode", y = NULL) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = .DG_title_size*0.8),
          axis.title      = element_text(size = .DG_axis_size),
          axis.text       = element_text(size = .DG_axis_size),
          legend.title    = element_text(size = .DG_legend_size),
          legend.text     = element_text(size = .DG_legend_size)
        )
     
      filter_comp <- gridExtra::arrangeGrob(
        plot_count_filter, plot_CR_filter,
        nrow = 1,
        top = grid::textGrob(
          "Comparison of blank-based filters for genes",
          gp = grid::gpar(fontsize = .DG_title_size)
          )
        )
     
      gridExtra::grid.arrange(
        grid::textGrob("PPV", rot = 90, gp = grid::gpar(fontsize = .DG_axis_size)),
        filter_comp,
        ncol = 2,
        widths = grid::unit.c(grid::unit(1, "lines"), grid::unit(1, "null"))
      )
     
    }
   
    return(results)
   
  }

# Plot run QC (PPV cutoff)
.DG_plot_PPV_cutoff <- function(
    PPV_resamples,
    CI_range = 0.95,
    min_PPV  = 0.8
  ) {
   
    CI_range_low  <- (1 - CI_range)/2
    CI_range_high <- CI_range + CI_range_low
   
    # Compute column medians
    lower_q      <- apply(PPV_resamples, 2, function(x) quantile(x, CI_range_low, na.rm = TRUE))
    # Order columns by decreasing lower quantile
    ordered_cols <- names(sort(lower_q, decreasing = TRUE))
   
    # Reorder columns in PPV_resamples
    PPV_resamples_ordered        <- PPV_resamples[ , ordered_cols]
    col_range                    <- sample(1:ncol(PPV_resamples_ordered), size = ncol(PPV_resamples_ordered), replace = FALSE)
    col_range                    <- sort(col_range, decreasing = FALSE)
    PPV_resamples_long           <- stack(PPV_resamples_ordered)
    colnames(PPV_resamples_long) <- c("value", "variable")
   
    # Summarize median and 95% confidence interval
   
    summary_df <- PPV_resamples_long |>
      dplyr::group_by(variable) |>
      dplyr::summarize(
        median = median(value, na.rm = TRUE),
        lower  = quantile(value, CI_range_low, na.rm = TRUE),
        upper  = quantile(value, CI_range_high, na.rm = TRUE)
      )
    summary_df$above_cutoff                             <- rep("bad", nrow(summary_df))
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
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        axis.title         = element_text(size = 18),               # Axis titles
        axis.text          = element_text(size = 12),               # Axis tick labels
        plot.title         = element_text(size = 24, hjust = 0.5),  # Title
        legend.text        = element_text(size = 12),
        legend.title       = element_text(size = 14)) +
      scale_color_manual(values = c("good" = "blue", "bad" = "black")) +
      guides(color = "none") +
      labs(
        y     = paste0("Estimated PPV: Median with ", CI_range, " CI"),
        x     = paste0("Gene barcodes by lower ", CI_range_low, " quantile of PPV confidence interval"),
        title = "Estimated Precision (PPV) by Gene Barcode")
    # ^ Save as 1600 x 925
   
    return(
      list(PPV_plot = PPV_plot, summary_df = summary_df)
    )
   
  }

# Plot simulation QC
.DG_plot_counts_sorted <- function(
    counts_sorted,
    obs_counts_sorted,
    N_genes
  ) {
   
    bc_type_colors <- c(
      "Gene (sim)" = .DG_retained_color,     "Blank (sim)" = .DG_blanks_color,
      "Gene (run)" = .DG_retained_color_obs, "Blank (run)" = .DG_blanks_color_obs
    )
   
    bc_type                                             <- rep("Gene (sim)", length(counts_sorted))
    bc_type[seq_along(counts_sorted) > N_genes]         <- "Blank (sim)"
    bc_type_obs                                         <- rep("Gene (run)", length(obs_counts_sorted))
    bc_type_obs[seq_along(obs_counts_sorted) > N_genes] <- "Blank (run)"
   
    df_counts_sorted <- data.frame(
      index   = c(seq_along(counts_sorted), seq_along(obs_counts_sorted)),
      Count   = c(counts_sorted, obs_counts_sorted),
      bc_type = c(bc_type, bc_type_obs),
      size    = c(rep(3, length(counts_sorted)), rep(1.5, length(obs_counts_sorted)))
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
        color      = "black",
        linetype   = "dashed",
        linewidth  = 0.5
      ) +
      labs(title = "Simulated vs Run Counts", x = "Barcodes sorted by run count", y = "Spot count", color = "Count type") +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title      = element_text(hjust = 0.5, size = .DG_title_size),
        axis.title      = element_text(size = .DG_axis_size),
        axis.text.y     = element_text(size = .DG_axis_size),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.title    = element_text(size = .DG_legend_size),
        legend.text     = element_text(size = .DG_legend_size)
      )
   
    return(counts_sorted_plot)
   
  }

# Plot simulation step information
.DG_plot_MCMC_steps <- function(
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

#' Run misread QC analysis using the legacy dichotomized-Gaussian (DG) spot simulation method
#'
#' \code{misread.qc.DG.legacy} is a verbatim-ported, historical alternative to \code{misread.qc}: rather
#' than fitting an analytic conditional-probability model, it estimates bit-flip rates/correlations (and
#' the resulting PPV) by directly simulating "ground truth" RNA hybridization spots and noisy MERFISH
#' barcode decoding, then adjusting the simulated flip rates and bit-noise correlations with a Markov
#' Chain Monte Carlo with Simulated Annealing (MCMCSA) search until the simulated barcode-count
#' distribution matches the observed one (summarized by \code{STdata}). This is the method used to
#' produce this project's earlier (pre-\code{misread.qc}) results; it is kept here, fully separate from
#' \code{misread.qc}/\code{sim.benchmark}, only so those results can in principle be reproduced, and as a
#' reference for \code{misread.qc.DG} (the C++-backed, actively-developed DG method) -- it is not the
#' package's primary (recommended) method, and is much slower than the analytic model for comparable
#' codebooks.
#'
#' Unlike \code{misread.qc}, \code{maxHam} here is not derived from the codebook; the original method
#' hard-coded a maximum correctable Hamming distance of 4 (an "extended" Hamming code assumption), which is
#' kept as the default here.
#'
#' @param STdata Numeric matrix with rows as barcodes, columns labeled "rates", "variance", "counts". Must have barcode names (e.g., gene or protein species) as row names, with blank species having "Blank" in their name.
#' @param codebook Numeric matrix with barcodes as row names and bits as columns. All entries should be 1 or 0. All row names from \code{STdata} must be included in the row names for \code{codebook}.
#' @param maxHam Maximum correctable Hamming distance used for error correction, default 4 (matches the original method's hard-coded value, an "extended" Hamming code assumption).
#' @param resamples Number of steps to take in the post-burn-in random walk (these are the steps summarized in the returned results).
#' @param burnin Number of burn-in steps before resampling begins, default \code{resamples/2}.
#' @param FR_search Number of initial steps used to search for a reasonable overall (mean) flip rate before burn-in, default \code{resamples/2}.
#' @param max_step_size Size of the largest per-step Gaussian proposal in the random walk, default 0.2.
#' @param min_step_size_divider Step size is scaled down linearly to \code{max_step_size/min_step_size_divider} over burn-in, default 10.
#' @param max_temperature Simulated-annealing temperature at the start of the walk, default 0.2.
#' @param min_temperature_scale Temperature is scaled down linearly to \code{max_temperature/(min_step_size_divider/min_temperature_scale)} over burn-in, default 0.5.
#' @param max_blank_weight Maximum weight given to blanks (vs. genes) in the fit objective during \code{FR_search}, default \code{N_genes/N_blanks + 1}.
#' @param min_blank_weight Minimum (final) weight given to blanks, default 2.
#' @param n_tracker_updates Number of progress updates to print over the full walk, default 100.
#' @param n_forks Number of process forks to use for simulating spots (i.e., parallel computation), default is 1. Must be 1 on Windows, can be higher on Mac and Linux.
#' @param N_cells Number of cells to simulate at each step, default 1000 (matches the original method's hard-coded default). Lowering this trades simulation accuracy for speed -- useful for quick/development runs.
#' @return A list giving:\itemize{
#'    \item \code{sim_summaries}: a data frame, one row per accepted MCMCSA step, of simulation-fit summary statistics (CR/PPV slope and MSE, stipulated flip rates/correlations, decoding rates, etc.).
#'    \item \code{PPV}: a data frame, one row per accepted MCMCSA step and one column per gene barcode, of the simulated PPV for that barcode at that step.
#'    \item \code{PPV_summary}: a data frame summarizing \code{PPV} (over the \code{resamples} rows following burn-in) per gene barcode: median, lower/upper 95\% CI, and whether the lower CI clears 0.8.
#'    \item \code{Counts_sorted_sim}/\code{Counts_sorted_obs}: data frames, one row per step, of the simulated and observed barcode counts (sorted consistently for comparison).
#'    \item \code{PPV_plot}: a ggplot object visualizing \code{PPV_summary}.
#'    \item \code{sorted_count_plot}: a ggplot object comparing simulated vs. observed sorted barcode counts over the \code{resamples} rows following burn-in.
#' }
#' @export
misread.qc.DG.legacy <- function(
    STdata,
    codebook,
    maxHam                = 4,     # Maximum correctable Hamming distance (original method hard-coded this to 4)
    resamples             = 1000,  # Number of steps to take in random walk
    burnin                = NULL,
    FR_search             = NULL,
    max_step_size         = 0.2,   # Size of each step in the random walk
    min_step_size_divider = 10,    # Scale step-size down linearly to max_step_size/min_step_size_divider over burnin
    max_temperature       = 0.2,   # Temperature for the MCMC simulation
    min_temperature_scale = 0.5,   # Scale temperature down linearly to max_temperature/(min_step_size_divider/min_temperature_scale) over burnin
    max_blank_weight      = NULL,
    min_blank_weight      = NULL,
    n_tracker_updates     = 100,
    n_forks               = 1,
    N_cells               = 1e3    # Number of cells to simulate per step (originally hard-coded via .DG_run_simulation()'s own default; exposed here, same default, so it can be lowered for quick/dev runs)
  ) {

    cat("\n\nEstimating FPQC metrics with Markov chain Monte Carlo simulated annealing (legacy DG method)")

    # Confirm forking is possible and check number of cores
    n_forks <- .DG_check_forks(n_forks)

    # Prep codebook and STdata, and split into the genes/blanks summary-stats data frames the original
    # method expects (mirrors the blank-mask/alignment logic the historical misreadQC() R wrapper used,
    # kept self-contained here rather than shared with main.R's check_data())
    species_names <- rownames(STdata)
    if (is.null(species_names))                      stop("STdata must have row names as species names")
    if (is.null(rownames(codebook)))                  stop("codebook must have row names as species names")
    if (!all(species_names %in% rownames(codebook)))  stop("All species in STdata must be present in codebook")
    blank_mask <- grepl("Blank", species_names, ignore.case = FALSE)
    if (sum(blank_mask) == 0)                         stop("No blanks found in STdata row names, make sure blank species have 'Blank' in their names")
    if (sum(blank_mask) == length(species_names))     stop("All species are blanks. Make sure non-blank species do not have 'Blank' in their names")
    summary_stats_genes  <- STdata[!blank_mask, , drop = FALSE]
    summary_stats_blanks <- STdata[ blank_mask, , drop = FALSE]
    # ... reorder codebook so gene rows precede blank rows, in the same order as the split summary stats
    #     above: .DG_run_simulation() matches summary-stats rows to codebook rows positionally (by order
    #     of appearance), mirroring how the original script's codebook/summary-stats CSVs were pre-aligned
    codebook <- codebook[c(rownames(summary_stats_genes), rownames(summary_stats_blanks)), , drop = FALSE]

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
    sim_results <- .DG_run_simulation(
      codebook = codebook,
      summary_stats_genes = summary_stats_genes,
      summary_stats_blanks = summary_stats_blanks,
      flip_rate_10 = flip_rate_10_current,
      flip_rate_01 = flip_rate_01_current,
      maxHam = maxHam,
      bit_lum_noise_correlations = lum_noise_correlations_current,
      blank_weight = max_blank_weight,
      maintain_gene_identity = FALSE,
      return_flip_rates = FALSE,
      n_forks = n_forks,
      N_cells = N_cells,
      verbose = FALSE
    )

    # Extract sim_mse of the initial simulation
    sim_mse_current <- sim_results$sim_summary["sim_mse"]
    cat("\nInitial mse:", round(sim_mse_current, 3))

    cat("\n\nRunning Markov chain steps: \n\nRunning initial steps ...")
    first_step <- TRUE
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
      sim_results <- .DG_run_simulation(
        codebook = codebook,
        summary_stats_genes = summary_stats_genes,
        summary_stats_blanks = summary_stats_blanks,
        flip_rate_10 = flip_rate_10_next,
        flip_rate_01 = flip_rate_01_next,
        maxHam = maxHam,
        bit_lum_noise_correlations = lum_noise_correlations_next,
        blank_weight = blank_weight_schedule[step + 1],
        maintain_gene_identity = maintain_gene_identity,
        return_flip_rates = return_flip_rates,
        n_forks = n_forks,
        N_cells = N_cells,
        verbose = FALSE
      )

      # Extract log-likelihood of the next simulation
      sim_mse_next <- sim_results$sim_summary["sim_mse"]

      # Use simulated annealing to decide whether to accept or reject the proposed step

      # Calculate acceptance probability
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
    PPV_plot_info <- .DG_plot_PPV_cutoff(PPV_genes[resamples_rows, ])
    PPV_summary <- PPV_plot_info$summary_df
    PPV_plot <- PPV_plot_info$PPV_plot
    print(PPV_plot)

    # Make sorted count plot to evaluate simulation
    counts_sorted_sim2 <- Counts_sorted_sim[resamples_rows,]
    counts_sorted_sim_ <- apply(counts_sorted_sim2, 2, median)
    counts_sorted_obs_ <- Counts_sorted_obs[1,]
    counts_sorted_obs_ <- unlist(as.vector(counts_sorted_obs_))
    sorted_count_plot <- .DG_plot_counts_sorted(counts_sorted_sim_, counts_sorted_obs_, N_genes)
    print(sorted_count_plot)

    # Plot steps to evaluate the resampling
    .DG_plot_MCMC_steps(sim_summaries, burnin, FR_search)

    # Bundle results
    results <- list(
      sim_summaries = sim_summaries,
      PPV = PPV_genes,
      PPV_summary = PPV_summary,
      Counts_sorted_sim = Counts_sorted_sim,
      Counts_sorted_obs = Counts_sorted_obs,
      PPV_plot = PPV_plot,
      sorted_count_plot = sorted_count_plot
    )

    return(results)

  }
