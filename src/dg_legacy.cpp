
// #########################################################################################################
// Verbatim port of paper_calling_rate_sim/MERFISH_calling_rate_sim/fast_functions.cpp, the C++ half of the
// original "dichotomized Gaussian" (DG) spot-simulation method used to produce the paper's results. Kept
// in its own file, with no changes to any function body, so those results can in principle be reproduced.
// Paired with R/dg_legacy.R (the R half of the same port). Fully separate from src/main.cpp's analytic
// model -- not wired together, not sharing any helper functions.
// #########################################################################################################

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector find_HD_by_row(
    const NumericMatrix& codebook,
    const NumericVector& spot
  ) {
    int n = codebook.nrow();
    NumericVector HD = NumericVector(n);
    for (int i = 0; i < n; i++) {
      HD[i] = sum(abs(codebook(i, _) - spot));
    }
    return HD;
  }

// [[Rcpp::export]]
List find_spot_labels_fast(
    const IntegerMatrix& codebook,
    const IntegerMatrix& spots,
    const int& max_correctable_Hamming_distance
  ) {
    int n_spots = spots.nrow();
    int n_codes = codebook.nrow();
    int m = codebook.ncol(); // dimension of each spot/codeword
    IntegerVector spot_labels(n_spots, NA_INTEGER);
    IntegerVector spot_labels_uncorrected(n_spots, NA_INTEGER);

    for (int i = 0; i < n_spots; i++) {
      int min_HD = m + 1;  // max possible HD is m
      int min_idx = -1;
      bool tie = false;

      for (int j = 0; j < n_codes; j++) {
        int h = 0;
        for (int k = 0; k < m; k++) {
          h += (static_cast<int>(spots(i, k)) != static_cast<int>(codebook(j, k)));
          if (h > min_HD) break;
        }
        if (h < min_HD) {
          min_HD = h;
          min_idx = j;
          tie = false;
        } else if (h == min_HD) {
          tie = true;
        }
        if (min_HD == 0) {
          break;
        }
      }

      if (min_HD <= max_correctable_Hamming_distance && !tie) {
        spot_labels[i] = min_idx + 1; // 1-indexed for R
        if (min_HD == 0) {
          spot_labels_uncorrected[i] = min_idx + 1; // 1-indexed for R
        }
      }

    }

    List spot_labels_list = List::create(
      Named("spot_labels") = spot_labels,
      Named("spot_labels_uncorrected") = spot_labels_uncorrected
    );

    return spot_labels_list;
  }

// [[Rcpp::export]]
NumericMatrix generate_spots_ci(
    IntegerMatrix bc_counts,
    NumericMatrix codebook_gen_mat,
    NumericMatrix noise,
    NumericMatrix corr_mat,
    int ci,
    double threshold,
    bool decode_and_label
  ) {
    int N_barcodes = bc_counts.ncol();
    int N_bits = codebook_gen_mat.ncol();

    int total_spots = 0;
    for (int j = 0; j < N_barcodes; j++) {
      total_spots += bc_counts(ci, j);
    }

    NumericMatrix spots_ci(total_spots, N_bits);
    std::fill(spots_ci.begin(), spots_ci.end(), NA_REAL);

    int end_row = 0;

    // Rcpp::Environment::namespace_env() loads MASS's namespace on demand (via R's asNamespace()),
    // unlike Environment("package:MASS") which requires MASS to already be attached to the search
    // path -- a package's Imports are loaded, not attached, so the original form would fail here.
    Environment MASS = Environment::namespace_env("MASS");
    Function mvrnorm = MASS["mvrnorm"];

    for (int bci = 0; bci < N_barcodes; bci++) {
      int count = bc_counts(ci, bci);
      end_row += count;
      int start_row = end_row - count;

      if (count >= 1) {
        // Get this barcode's mean vector
        NumericVector this_bc(N_bits);
        for (int k = 0; k < N_bits; k++) {
          this_bc[k] = codebook_gen_mat(bci, k);
        }

        // Build noise diagonal matrix
        NumericMatrix noise_matrix(N_bits, N_bits);
        for (int i = 0; i < N_bits; i++) {
          for (int j = 0; j < N_bits; j++) {
            if (i == j) {
              noise_matrix(i, j) = noise(bci, i);
            } else {
              noise_matrix(i, j) = 0.0;
            }
          }
        }

        // Sigma_noised = noise_matrix %*% corr_mat %*% noise_matrix
        NumericMatrix temp_sigma(N_bits, N_bits);
        for (int i = 0; i < N_bits; i++) {
          for (int j = 0; j < N_bits; j++) {
            double sum = 0.0;
            for (int k = 0; k < N_bits; k++) {
              sum += noise_matrix(i, k) * corr_mat(k, j);
            }
            temp_sigma(i, j) = sum;
          }
        }

        NumericMatrix Sigma_noised(N_bits, N_bits);
        for (int i = 0; i < N_bits; i++) {
          for (int j = 0; j < N_bits; j++) {
            double sum = 0.0;
            for (int k = 0; k < N_bits; k++) {
              sum += temp_sigma(i, k) * noise_matrix(k, j);
            }
            Sigma_noised(i, j) = sum;
          }
        }

        // Sample from multivariate normal using MASS::mvrnorm
        NumericMatrix spot_seeds(count, N_bits);
        if (count > 1) {
          spot_seeds = mvrnorm(count, this_bc, Sigma_noised);
        } else {
          spot_seeds.row(0) = Rcpp::as<NumericVector>(mvrnorm(count, this_bc, Sigma_noised));
        }

        for (int row = 0; row < count; row++) {
          for (int col = 0; col < N_bits; col++) {
            double val = spot_seeds(row, col);
            if (decode_and_label) {
              spots_ci(start_row + row, col) = (val > threshold) ? 1.0 : 0.0;
            } else {
              spots_ci(start_row + row, col) = val;
            }
          }
        }
      }
    }

    return spots_ci;
  }

// [[Rcpp::export]]
IntegerVector find_nearby_barcodes_fast(
  const int& bci,
  const IntegerMatrix& codebook,
  const int& max_correctable_Hamming_distance
  ) {
    int n_barcodes = codebook.nrow();
    int n_bits = codebook.ncol();
    IntegerVector barcode = codebook(bci, _);
    IntegerVector nearby_barcodes;
    for (int i = 0; i < n_barcodes; i++) {
      IntegerVector codebook_i = codebook(i, _);
      int HD = 0;
      for (int k = 0; k < n_bits; k++) {
        HD += (static_cast<int>(codebook_i[k]) != static_cast<int>(barcode[k]));
        if (HD > max_correctable_Hamming_distance) {
          break; // No need to continue if HD exceeds the limit
        }
      }
      if (HD <= max_correctable_Hamming_distance) {
        nearby_barcodes.push_back(i + 1); // R is 1-indexed
      }
    }
    return(nearby_barcodes);
  }

// Function to compute PPV for each barcode
// ... want to answer the question: for each barcode bc, if a spot is decoded as bc, what is the probability this is correct?
// [[Rcpp::export]]
NumericVector compute_PPV_fast(
    const List& barcodes_by_cell_true_list,
    const List& spots_bc_labels_list,
    const int& N_barcodes
  ) {

    // Initialize vectors to hold counts
    int N_cells = barcodes_by_cell_true_list.size();
    IntegerVector decoded_counts(N_barcodes, 0);
    IntegerVector successful_decoded_counts(N_barcodes, 0);

    for (int ci = 0; ci < N_cells; ci++) {
      // Get decoded counts and true barcode IDs for this cell
      IntegerVector true_barcodes_ids = barcodes_by_cell_true_list[ci];
      IntegerVector decoded_barcode = spots_bc_labels_list[ci];
      LogicalVector successfully_decoded = !is_na(decoded_barcode);
      int n_spots = true_barcodes_ids.size();
      for (int si = 0; si < n_spots; si++) {
        if (successfully_decoded[si]) {
          // For this spot, get its true barcode
          int true_bc = true_barcodes_ids[si] - 1; // R is 1-indexed, convert to 0-indexed
          // Get the decoded barcode
          int this_read = decoded_barcode[si] - 1; // R is 1-indexed, convert to 0-indexed
          // Increment the counts
          decoded_counts[this_read]++;
          if (this_read == true_bc) {
            successful_decoded_counts[this_read]++;
          }
        }
      }
    }

    NumericVector PPV_true = as<NumericVector>(successful_decoded_counts) / as<NumericVector>(decoded_counts);
    return PPV_true;

  }
