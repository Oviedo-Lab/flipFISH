
// Rcpp
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <cmath>
#include <unistd.h>     // fork, pipe, read, write, close
#include <sys/wait.h>   // waitpid
using namespace Rcpp;
using namespace Eigen;

// Hamming distance
// int dist = __builtin_popcountll(a ^ b);

/*
 * *********************************************************************************************************************
 * Basic data structures
 */

struct Codebook {
  int N_bits;
  std::vector<uint64_t> barcodes;
  std::vector<std::string> species;
  std::vector<int> blanks;
  std::vector<int> genes; 
};

struct EvalHist {
  std::vector<double> msle;
  std::vector<std::vector<double>> fr;
  std::vector<std::vector<double>> etc;
  std::vector<std::vector<double>> ecc; 
  std::vector<std::vector<double>> erc;
  std::vector<std::vector<double>> CR;
  std::vector<std::vector<double>> PPV;
  std::vector<std::vector<double>> est_true_bc_counts;
};

struct ST_data {
  // ST data
  std::vector<double> bc_rates;           // Expected count, per cell
  std::vector<double> bc_variance;        // Expected variance in count, among cells
  std::vector<int> bc_counts;             // Vector of length N_barcodes, giving total counts for each barcode across all cells
  // Codebook
  Codebook cb;
  int max_correctable_Hamming_distance;
  std::unordered_map<uint64_t, int> correction_table;
  std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted;
  // Parameter estimation
  EvalHist eval_history;
  int n_forks;
  int report_freq;
};

struct FlipRates {
  std::vector<double> rate10; // P(1 -> 0) for each bit, so length = N_bits
  std::vector<double> rate01; // P(0 -> 1) for each bit, so length = N_bits
  std::vector<double> log_rate10; // Pre-computed vector of log(rate10) values, used in TR calculation
  std::vector<double> log_rate01; // Pre-computed vector of log(rate01) values, used in TR calculation
  std::vector<double> log_rate11; // Pre-computed vector of 1-log(rate10) values, used in TR calculation
  std::vector<double> log_rate00; // Pre-computed vector of 1-log(rate01) values, used in TR calculation
  std::vector<double> corr1; // Strict lower triangle of luminance noise correlations when bit i is 1 
  std::vector<double> corr0; // ... when bit i is 0
  std::vector<std::vector<double>> log_inv_corr; // Pre-computed vector of summed j<i values log(1 - corr(i,j)), used in TR calculation
}; 

/*
 * *********************************************************************************************************************
 * Helper functions, mathematical
 */

MatrixXd make_corr_matrix(
    const std::vector<double>& corr_flat, 
    int N_bits
  ) {
    MatrixXd L = MatrixXd::Zero(N_bits, N_bits);
    int idx = 0;
    for (int i = 0; i < N_bits; ++i) {
        // positive diagonal
        L(i,i) = std::exp(corr_flat[idx++]);
        for (int j = 0; j < i; ++j) {L(i,j) = corr_flat[idx++];}
    }
    MatrixXd S = L * L.transpose();
    VectorXd d = S.diagonal().array().sqrt();
    MatrixXd corr(N_bits, N_bits);
    for (int i = 0; i < N_bits; ++i) {
      for (int j = 0; j < N_bits; ++j) {corr(i,j) = S(i,j) / (d(i) * d(j));}
    }
    return corr;
  }

double compute_msle(
    const std::vector<int>& obs_counts,
    const std::vector<double>& pred_counts,           // predicted corrected counts
    const std::vector<double>& weights
  ) {
    int N_barcodes = obs_counts.size();
    if (pred_counts.size() != N_barcodes || weights.size() != N_barcodes) {Rcpp::stop("obs_counts, pred_counts, and weights must be the same length.");}
    double msle = 0.0;
    for (int b = 0; b < N_barcodes; ++b) {
      double sle = std::log(pred_counts[b] + 1.0) - std::log((double)obs_counts[b] + 1.0);
      msle += sle * sle * weights[b];
    }
    msle /= (double)N_barcodes;
    return msle;
  }

// Formula to calculate gamme dispersion factor from mean and variance of counts
// ... taken from the wispack code
double compute_gamma_dispersion(
    const double& count_mean, // mean of counts for context-species combination
    const double& count_var   // variance of counts for context-species combination
  ) {
    if (count_var > count_mean) {
      // Have: count_var = count_mean + count_mean^2 * gdis
      // ... count_var = count_mean * (1 + count_mean * gdis)
      // ... count_var / count_mean = 1 + count_mean * gdis
      // ... (count_var / count_mean) - 1 = count_mean * gdis
      return ((count_var / count_mean) - 1.0) / count_mean;
    } else {
      return 0.0; // no dispersion if variance is less than or equal to mean
    }
  }

// Integral of Poisson-Gamma distribution, from 0 to positive infinity
// ... taken from the wispack code
double poisson_gamma_integral(
    double y, // observed count value
    double r, // expected process rate drawn from the gamma distribution
    double v  // variance of the gamma distribution
  ) {
   
    // convert more intuitive rate and variance parameters into the standard shape-rate parameters 
    // ... for the gamma distribution
    double s = r * r / v; // Gamma distribution "shape" parameter
    double R = s / r;     // Gamma distribution "rate" parameter
    
    // Idea: This is an analytic solution to the integral of dPois(y, lambda) * dGamma(lambda, r, v) from 0 to positive infinity. 
    //        ... The solution takes the form of a ratio num/denom which is subject to overflow/underflow, and so instead of 
    //            computing this ratio directly, we compute the log of the numerator and denominator separately, and then exponentiate the difference.
    double log_num = s * std::log(R) + std::lgamma(y + s); 
    double log_denom = std::lgamma(y + 1.0) + std::lgamma(s) + (y + s) * std::log(R + 1.0);
    double integral = std::exp(log_num - log_denom);
   
    // note: this return value is *not* the log of the density! It's the density itself. 
    return integral;
    
  }

double compute_nll(
    const ST_data& STdata,
    const std::vector<double>& pred_counts           // predicted corrected counts
  ) {
    // Compute negative log likelihood of bc_counts, given corrected_counts as predicted counts
    double log_lik = 0.0;
    double N_cells = (double)STdata.bc_counts[0] / STdata.bc_rates[0]; // Estimate number of cells from observed counts and rates for one of the barcodes
    for (int i = 0; i < pred_counts.size(); ++i) {
      double pred_rate = pred_counts[i] / N_cells;
      if (pred_rate == 0.0) {pred_rate = 1e-12;} // Avoid zero predicted rates, which cause issues for the likelihood calculation
      double gamma_variance = pred_rate + compute_gamma_dispersion(STdata.bc_rates[i], STdata.bc_variance[i]) * pred_rate * pred_rate;
      log_lik += std::log(poisson_gamma_integral(STdata.bc_rates[i], pred_rate, gamma_variance));
    }
    return -log_lik;
  }

/*
 * *********************************************************************************************************************
 * Helper functions, data manipulation
 */

std::vector<int> grep_idx(
    Rcpp::CharacterVector x,
    std::string s,
    bool neg = false
  ) {
    std::vector<int> idx;
    for (int i = 0; i < x.size(); ++i) {
      std::string xi = Rcpp::as<std::string>(x[i]);
      if (xi.find(s) != std::string::npos) {
        if (!neg) {
          idx.push_back(i);  
        }
      } else if (neg) {
        idx.push_back(i);
      }
    }
    if (idx.size() > 0) {
      return idx; 
    } else {
      return {-1};
    }
  } 

std::vector<double> compute_log_inv_corr(
    const FlipRates& fr,
    uint64_t bc,
    uint64_t flips,
    int N_bits
  ) {
    std::vector<double> log_inv_corr(N_bits, 0.0);
    std::vector<int> flip(N_bits, 0);
    // Compute log of the inverse of the correlation between bit flips
    for (int i = 0; i < N_bits; ++i) {
      int bit = (bc >> i) & 1ULL;
      flip[i] = (flips >> i) & 1ULL;
      for (int j = 0; j < i; ++j) {
        if (bit & flip[j]) {
          log_inv_corr[i] += std::log(1.0 - fr.corr1[i * (i - 1) / 2 + j]);
        } else if (flip[j]) {
          log_inv_corr[i] += std::log(1.0 - fr.corr0[i * (i - 1) / 2 + j]);
        }
      }
    }
    return log_inv_corr;
  }

FlipRates pack_fr(
    const std::vector<double>& rates_plus_corr,
    int N_bits
  ) {
    // Initialize flip-rate structure
    FlipRates fr;
    // Assign flip rates
    fr.rate10.assign(rates_plus_corr.begin(), rates_plus_corr.begin() + N_bits);
    fr.rate01.assign(rates_plus_corr.begin() + N_bits, rates_plus_corr.begin() + 2*N_bits);
    // Find number of correlations (strict lower triangle) ... corr(i,j) = corrXXXX(i * (i - 1) / 2 + j)
    int N_corrs = N_bits * (N_bits - 1) / 2;
    // Assign correlations
    fr.corr1.assign(rates_plus_corr.begin() + 2*N_bits, rates_plus_corr.begin() + 2*N_bits + N_corrs);
    fr.corr0.assign(rates_plus_corr.begin() + 2*N_bits + N_corrs, rates_plus_corr.end());
    // Precompute log of the rates and log of 1 - rates, which are used in the TR calculation
    fr.log_rate00 = std::vector<double>(N_bits, 0.0);
    fr.log_rate01 = std::vector<double>(N_bits, 0.0);
    fr.log_rate10 = std::vector<double>(N_bits, 0.0);
    fr.log_rate11 = std::vector<double>(N_bits, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      fr.log_rate10[i] = std::log(fr.rate10[i]);
      fr.log_rate01[i] = std::log(fr.rate01[i]);
      fr.log_rate00[i] = std::log(1.0 - fr.rate01[i]);
      fr.log_rate11[i] = std::log(1.0 - fr.rate10[i]);
    }
    return fr;
  }

std::unordered_map<int, std::vector<uint64_t>> invert_lookup_table(
    const std::unordered_map<uint64_t, int>& lut
  ) {
    std::unordered_map<int, std::vector<uint64_t>> inverted;
    for (const auto& kv : lut) {
      uint64_t barcode = kv.first;
      int label = kv.second;
      inverted[label].push_back(barcode);
    }
    return inverted;
  }

uint64_t pack(
    std::vector<int> bits
  ) {
    uint64_t packed = 0;
    for (int i = 0; i < bits.size(); ++i) {packed |= (uint64_t(bits[i]) << i);} 
    return packed;
  }

void generate_neighbors(
    uint64_t x,
    int n_bits,
    int dist,
    int start_bit,
    std::vector<uint64_t>& out
  ) {
    if (dist == 0) {
      out.push_back(x);
      return;
    }
    for (int b = start_bit; b < n_bits; ++b) {
      uint64_t flipped = x ^ (1ULL << b); // Flip bit b
      generate_neighbors(
        flipped,
        n_bits,
        dist - 1,
        b + 1,
        out
      );
    }
  }

std::vector<uint64_t> neighbors(
    uint64_t x,
    int n_bits,
    int hamming_dist
  ) {
    std::vector<uint64_t> out;
    generate_neighbors(
      x,
      n_bits,
      hamming_dist,
      0, // start bit
      out
    );
    return out;
  }

/*
 * *********************************************************************************************************************
 * Functions for spot decoding
 */

std::vector<int> unique_Hamming(
    const std::vector<uint64_t>& barcodes
  ) {
    int N_barcodes = barcodes.size();
    std::vector<int> hamming_distances;
    for (int i = 0; i < N_barcodes; ++i) {
      for (int j = i + 1; j < N_barcodes; ++j) {
        int dist = __builtin_popcountll(barcodes[i] ^ barcodes[j]);
        hamming_distances.push_back(dist);
      }
    }
    std::sort(hamming_distances.begin(), hamming_distances.end());
    auto last = std::unique(hamming_distances.begin(), hamming_distances.end());
    hamming_distances.erase(last, hamming_distances.end());
    return hamming_distances;
  }

Codebook pack_codebook(
    const IntegerMatrix& codebook,
    bool verbose = true
  ) {
    // Extract basic info
    int N_bits = codebook.ncol();
    int N_barcodes = codebook.nrow();
    CharacterVector species = rownames(codebook);
    // Initialize new Codebook
    Codebook cb;
    cb.N_bits = N_bits;
    // Make index for blanks and genes
    cb.blanks = grep_idx(species, "Blank");
    cb.genes = grep_idx(species, "Blank", true);
    // Extract barcodes as packed bits
    cb.barcodes.reserve(N_barcodes);
    cb.species.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      std::vector<int> bits(N_bits);
      for (int j = 0; j < N_bits; ++j) {bits[j] = codebook(i, j);}
      cb.barcodes.push_back(pack(bits));
      cb.species.push_back(as<std::string>(species[i]));
    }
    // Find all unique Hamming distances between barcodes 
    std::vector<int> hamming_distances = unique_Hamming(cb.barcodes);
    // Print codebook stats
    if (verbose) {
      Rcpp::Rcout << "Codebook stats:\n";
      Rcpp::Rcout << "  N_bits: " << N_bits << "\n";
      Rcpp::Rcout << "  N_barcodes: " << N_barcodes << "\n";
      Rcpp::Rcout << "  N_blanks: " << cb.blanks.size() << "\n";
      Rcpp::Rcout << "  Unique Hamming distances between barcodes: ";
      for (int d : hamming_distances) {Rcpp::Rcout << d << " ";}
      Rcpp::Rcout << std::endl;
    }
    // Return packed codebook
    return cb;
  }

// [[Rcpp::export]]
IntegerVector unique_Hamming_cb(
    const IntegerMatrix& codebook
  ) {
    Codebook cb = pack_codebook(codebook, false);
    std::vector<int> hamming_distances = unique_Hamming(cb.barcodes);
    return wrap(hamming_distances);
  }

std::unordered_map<uint64_t, int> build_correction_table(
    const Codebook& cb,
    int max_correctable_Hamming_distance
  ) {
    std::unordered_map<uint64_t, int> correction_table;
    for (size_t i = 0; i < cb.barcodes.size(); ++i) {
      correction_table[cb.barcodes[i]] = i; // Exact match
      // Generate all barcodes within max_correctable_Hamming_distance
      for (int d = 1; d <= max_correctable_Hamming_distance; ++d) {
        std::vector<uint64_t> n = neighbors(cb.barcodes[i], cb.N_bits, d);
        for (int j = 0; j < n.size(); ++j) {
          if (correction_table.count(n[j])) {
            // If this neighbor is already mapped to a different barcode, we have a tie
            if (correction_table[n[j]] != i) {
              correction_table[n[j]] = -2; // Mark as ambiguous
            }
          } else {
            correction_table[n[j]] = i; // Map neighbor to original barcode index
          }
        }
      }
    }
    return correction_table;
  }

/*
 * *********************************************************************************************************************
 * Functions for data loading 
 */

ST_data load_STdata(
    NumericMatrix bc_count_data,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance
  ) {
    
    // Load codebook as packed integers
    Codebook cb = pack_codebook(codebook);
    
    // Get info from bc_count_data 
    int N_barcodes = bc_count_data.nrow();
    int rate_col = -1;
    int variance_col = -1; 
    int count_col = -1;
    CharacterVector data_colnames = colnames(bc_count_data);
    for (int i = 0; i < bc_count_data.ncol(); ++i) {
      if (data_colnames[i] == "rates") {rate_col = i;} 
      else if (data_colnames[i] == "variance") {variance_col = i;} 
      else if (data_colnames[i] == "counts") {count_col = i;} 
    }
    if (rate_col == -1 || variance_col == -1 || count_col == -1) {
      Rcpp::stop("bc_count_data must have columns named 'rates', 'variance', and 'counts'.");
    }
    
    // Extract rate, variance, and count data into vectors
    std::vector<double> bc_rates;
    std::vector<double> bc_variance;
    std::vector<int> bc_counts;
    bc_rates.reserve(N_barcodes);
    bc_variance.reserve(N_barcodes);
    bc_counts.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      bc_rates.push_back((double)bc_count_data(i, rate_col));
      bc_variance.push_back((double)bc_count_data(i, variance_col));
      bc_counts.push_back((int)bc_count_data(i, count_col));
    }
    
    // Build correction table 
    std::unordered_map<uint64_t, int> correction_table = build_correction_table(cb, max_correctable_Hamming_distance);
    Rcpp::Rcout << "Correction table built with max-correctable Hamming distance " << max_correctable_Hamming_distance << "." << std::endl;
    Rcpp::Rcout << "Correction table size: " << correction_table.size() << " correctable barcodes." << std::endl;
    
    // Invert correction table
    std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted = invert_lookup_table(correction_table);
    Rcpp::Rcout << "Correction table inverted." << std::endl;
    
    // Return parsed ST data
    EvalHist eval_history;
    return {
      bc_rates, bc_variance, bc_counts, 
      cb, 
      max_correctable_Hamming_distance,
      correction_table, correction_table_inverted,
      eval_history, 0, 0 // placeholders
    };
    
  }

/*
 * *********************************************************************************************************************
 * Functions to analytically compute expected count read from flip rates and true counts
 */

// Log of the rate at which the sequence of flips transform_flips can be expected to occur for a spot with true barcode bc, given bit-flip rates rate10 and rate01 and correlation corr between bit-flips
// STOPPED ... rates need to be normalized so they add to 1 once correlations are included ?!?!?!?!
double TR(
    const uint64_t bc,
    const uint64_t transform_flips,
    const FlipRates& fr
  ) {
    // Get number of bits
    int N_bits = fr.rate10.size();
    // Initialize variable to hold log of transformation rate
    double log_tr = 0.0;
    std::vector<double> log_inv_corr = compute_log_inv_corr(fr, bc, transform_flips, N_bits);
    for (int i = 0; i < N_bits; ++i) {
      // Extract bit i from barcode bc
      int bit = (bc >> i) & 1ULL;
      // Extract whether bit i is flipped from transform_flips
      int flip = (transform_flips >> i) & 1ULL;
      // Compute independent transformation rate for this bit 
      double log_tr_nocorr = fr.log_rate00[i];
      if (flip) {
        if (bit) {log_tr_nocorr = fr.log_rate10[i];} 
        else {log_tr_nocorr = fr.log_rate01[i];}
      } else if (bit) {
        log_tr_nocorr = fr.log_rate11[i];
      }
      // Add adjusted log transformation rate for this bit to total log transformation rate
      log_tr += log_tr_nocorr - log_inv_corr[i];
    }
    return std::exp(log_tr); 
  }

// Function to compute the expected count for a given barcode of interest (BOI)
std::tuple<double, double, double> expected_bc_count(
    const uint64_t BOI,                             // Barcode of interest (BOI) for which we want to compute expected count after correction
    const FlipRates& fr,                            // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>& bc_counts_true,         // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,     // Vector giving all possible true spot barcodes
    const std::vector<uint64_t>& corrected_to_BOI   // Vector of barcodes that would be corrected to barcode of interest (BOI)
  ) {
    int N_bits = fr.rate10.size();
    int N_barcodes = true_barcodes.size();
    int N_correctable = corrected_to_BOI.size();
    double count_corrected = 0.0;
    double count_read = 0.0;
    double count_true = 0.0;
    // For each possible barcode misread that would be corrected to barcode BOI ...
    for (int j = 0; j < N_barcodes; ++j) {
      if (bc_counts_true[j] > 0) {
        double tr = 0.0;
        // Populate transformation rates
        std::vector<double> TRj(N_correctable);
        for (int i = 0; i < N_correctable; ++i) {
          // Get bit-flips required to transform the true barcode into the misread barcode
          uint64_t flips_to_BOI = (true_barcodes[j] ^ corrected_to_BOI[i]) & ((1ULL << N_bits) - 1);
          // Find expected count of spot misreads corrected to barcode BOI, from true barcode j, and add to total misread count
          TRj[i] = TR(true_barcodes[j], flips_to_BOI, fr);
        }
        // Sum transformation rates
        double sum_TRj = std::accumulate(TRj.begin(), TRj.end(), 0.0);
        if (sum_TRj > 0.0) {
          for (int i = 0; i < N_correctable; ++i) {
            // Normalize TRj to sum to 1
            TRj[i] /= sum_TRj;
            // Add to the transformation rate
            tr += TRj[i];
            if (corrected_to_BOI[i] == BOI) {count_read += TRj[i] * (double)bc_counts_true[j];}
          }
        } else {
          Rcpp::stop("Error in computing transformation rates.");
        }
        count_corrected += tr * (double)bc_counts_true[j];
        if (true_barcodes[j] == BOI) {count_true += tr * (double)bc_counts_true[j];}
      }
    }
    return {count_read, count_corrected, count_true};
  }

// Estimate expected barcode counts after correction, as a function of flip rates and true barcode counts, for all barcodes
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> expected_bc_counts(
    const FlipRates& fr,                                                              // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>& bc_counts_true,                                           // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,                                       // Vector giving all possible true spot barcodes
    const std::unordered_map<int, std::vector<uint64_t>>& correction_table_inverted,  // Inverted correction table mapping each barcode index to vector of misread barcodes that would be corrected to it
    int n_forks
  ) {
    int N_barcodes = true_barcodes.size();
    std::vector<double> ecc(N_barcodes, 0.0); // expected corrected counts
    std::vector<double> erc(N_barcodes, 0.0); // expected read counts
    std::vector<double> etc(N_barcodes, 0.0); // expected true (i.e., correctly read) counts
    
    // Construct batches of barcodes to simulate in parallel
    std::vector<std::vector<int>> barcode_batches(n_forks);
    for (int batch = 0; batch < n_forks; ++batch) {
      int start = (batch * N_barcodes) / n_forks;
      int end   = ((batch + 1) * N_barcodes) / n_forks;
      for (int i = start; i < end; ++i) {barcode_batches[batch].push_back(i);}
    }
    if (n_forks == 1 && barcode_batches[0].size() != N_barcodes) {
      Rcpp::stop("Error in batching barcodes for simulation. Expected all barcodes to be in one batch when n_forks=1.");
    }
    
    if (n_forks > 1) {
      // Run in parallel with forking
      
      // Pipes for inter-process communication
      std::vector<int> pids(n_forks);
      std::vector<std::array<int, 2>> pipes(n_forks); 
      
      // Initialize pipes 
      for (int i = 0; i < n_forks; ++i) {
        pipe(pipes[i].data());
      }
      
      // fork processes
      for (int i = 0; i < n_forks; i++) {
        
        pid_t pid = fork();
        int N_barcodes_batch = barcode_batches[i].size();
        
        if (pid == 0) { // child process
          
          // Close unrelated pipe ends
          for (int j = 0; j < n_forks; ++j) {
            if (j == i) {
              close(pipes[j][0]); // keep write end
            } else { 
              close(pipes[j][0]);
              close(pipes[j][1]);
            } 
          }
          std::vector<double> erctc_child(3*N_barcodes_batch, 0.0);
          
          for (int b = 0; b < N_barcodes_batch; ++b) {
            std::tuple<double, double, double> erctc = expected_bc_count(
              true_barcodes[barcode_batches[i][b]],
              fr, bc_counts_true, true_barcodes, 
              correction_table_inverted.at(barcode_batches[i][b])
            );
            erctc_child[b] = std::get<0>(erctc); // expected read count
            erctc_child[b + N_barcodes_batch] = std::get<1>(erctc); // expected corrected count
            erctc_child[b + 2*N_barcodes_batch] = std::get<2>(erctc); // expected true count
          } 
          
          // Send result 
          const char* buffer = reinterpret_cast<const char*>(erctc_child.data());
          size_t nbytes = sizeof(double) * erctc_child.size();
          size_t total_written = 0;
          while (total_written < nbytes) {
            ssize_t n_written = write(
              pipes[i][1],
              buffer + total_written,
              nbytes - total_written 
            );
            if (n_written <= 0) {
              if (errno == EINTR) continue; // Retry if interrupted
              close(pipes[i][1]);
              _exit(1);
            } 
            total_written += static_cast<size_t>(n_written);
          } 
          
          close(pipes[i][1]);    // Close write end
          _exit(0);              // Exit child process
          
        } else if (pid > 0) { // parent process 
          pids[i] = pid;      // Grab child pid
          close(pipes[i][1]); // Close write end
        } else { 
          Rcpp::stop("Fork failed!");
        } 
        
      } 
      
      // Fetch results from pipes
      for (int i = 0; i < n_forks; i++) {
        int N_barcodes_batch = barcode_batches[i].size();
        std::vector<double> erctc_child(3*N_barcodes_batch, 0.0);
        
        // Read the row from the pipe into the buffer
        char* buffer = reinterpret_cast<char*>(erctc_child.data());
        size_t nbytes = sizeof(double) * erctc_child.size();
        size_t total_read = 0;
        while (total_read < nbytes) {
          ssize_t n_read = read(
            pipes[i][0],
            buffer + total_read,
            nbytes - total_read
          );
          if (n_read == 0) {
            close(pipes[i][0]);
            int status;
            waitpid(pids[i], &status, WNOHANG);
            Rcpp::stop(
              "\nUnexpected EOF while reading pipe " + std::to_string(i) + " of " + std::to_string(n_forks) + ". " +
              "\nCheck for quiet core crash from a child-process error." +
              "\nRead " + std::to_string(total_read) + " of " + std::to_string(nbytes) + " bytes." + 
              "\nBatch size: " + std::to_string(barcode_batches[i].size()) + " barcodes, from " + std::to_string(barcode_batches[i][0]) + " to " + std::to_string(barcode_batches[i].back()) + "."
            );
          }
          if (n_read < 0) {
            if (errno == EINTR) continue; // Retry if interrupted
            close(pipes[i][0]);
            Rcpp::stop("Pipe read failed");
          }
          total_read += static_cast<size_t>(n_read);
        }
        
        for (int b = 0; b < N_barcodes_batch; ++b) {
          erc[barcode_batches[i][b]] += erctc_child[b];
          ecc[barcode_batches[i][b]] += erctc_child[b + N_barcodes_batch];
          etc[barcode_batches[i][b]] += erctc_child[b + 2*N_barcodes_batch];
          }
        close(pipes[i][0]);           // Close read end
        int status;
        waitpid(pids[i], &status, 0);
        if (WIFSIGNALED(status)) {Rcpp::Rcout << "Child killed by signal " << WTERMSIG(status) << std::endl;}
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {Rcpp::Rcout << "Child exited with code " << WEXITSTATUS(status) << std::endl;}
      }
      
    } else {
      // Run in serial
      for (int b : barcode_batches[0]) {
        std::tuple<double, double, double> erctc = expected_bc_count(
          true_barcodes[barcode_batches[0][b]],
          fr, bc_counts_true, true_barcodes, 
          correction_table_inverted.at(barcode_batches[0][b])
        );
        erc[b] = std::get<0>(erctc); // expected read count
        ecc[b + N_barcodes] = std::get<1>(erctc); // expected corrected count
        etc[b + 2*N_barcodes] = std::get<2>(erctc); // expected true counts
      }
    }
    
    return {erc, ecc, etc};
  }

/*
 * *********************************************************************************************************************
 * Functions to estimate flip rates, PPV, and CR
 */

// Estimate true barcode counts from observed counts plus flip rates
// ... merely intended to prevent msle from rising due to lack of spots
std::vector<int> est_bc_counts_true(
    const FlipRates& fr,
    void* data
  ) {
    // Grab data
    auto* d = static_cast<ST_data*>(data);
    const auto& true_barcodes = d->cb.barcodes;
    int N_barcodes = true_barcodes.size();
    int N_bits = d->cb.N_bits;
    
    // Find expected decoding rate
    double mean_Hamming_weight = 0.0; 
    for (int i = 0; i < N_barcodes; ++i) {mean_Hamming_weight += (double)__builtin_popcountll(true_barcodes[i]);}
    mean_Hamming_weight /= (double)N_barcodes;
    // ... find expected flip rate
    double mean_rate10 = std::accumulate(fr.rate10.begin(), fr.rate10.end(), 0.0);
    mean_rate10 /= (double)N_bits;
    double mean_rate01 = std::accumulate(fr.rate01.begin(), fr.rate01.end(), 0.0);
    mean_rate01 /= (double)N_bits;
    double expected_flip_rate = (mean_rate10*mean_Hamming_weight + mean_rate01*((double)N_bits - mean_Hamming_weight)) / (double)N_bits;
    // ... find degree of freedom
    double mean_corr1 = std::accumulate(fr.corr1.begin(), fr.corr1.end(), 0.0);
    mean_corr1 /= (double)fr.corr1.size();
    double mean_corr0 = std::accumulate(fr.corr0.begin(), fr.corr0.end(), 0.0);
    mean_corr0 /= (double)fr.corr0.size();
    double DoF = (1.0 - mean_corr1)*mean_Hamming_weight + (1.0 - mean_corr0)*((double)N_bits - mean_Hamming_weight);
    // ... find expected decoding rate
    double expected_decoding_rate = R::ppois(
      d->max_correctable_Hamming_distance, // Max number of bits that can flip and still be decoded
      expected_flip_rate * DoF, // Expected number of bit flips per spot, given the expected flip rate and number of bits
      true, false
      );
    
    // Adjust bc_counts for expected_decoding_rate
    std::vector<int> bc_counts_true(N_barcodes, 0);
    for (int i = 0; i < N_barcodes; ++i) {
      bc_counts_true[i] = (int)std::round((double)d->bc_counts[i] / expected_decoding_rate);
    }
   
    // Set blanks to zero
    for (int idx : d->cb.blanks) {
      bc_counts_true[idx] = 0;
    }
    
    return bc_counts_true; 
  }

// Make initial estimate of flip rates and bit-flip correlations
std::tuple<
    std::vector<double>, 
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
  > scale_initial_fr(
    const ST_data& STdata
  ) {
    Rcpp::Rcout << "\nEstimating scale of initial flip rates and correlations:" << std::endl;
    // Get codebook info
    int N_genes = STdata.cb.genes.size();
    int N_blanks = STdata.cb.blanks.size();
    int N_bits = STdata.cb.N_bits;
    int corr_free = N_bits * (N_bits - 1) / 2;
    // Initialize data structures to hold output
    std::vector<double> initial_rate10(N_bits, 0.0);
    std::vector<double> initial_rate01(N_bits, 0.0);
    std::vector<double> initial_corr1(corr_free, 0.0);
    std::vector<double> initial_corr0(corr_free, 0.0);
    // For each gene ...
    for (int gi = 0; gi < N_genes; ++gi) {
      // ... extract barcode bits
      int g = STdata.cb.genes[gi];
      std::vector<int> bit(N_bits);
      for (int i = 0; i < N_bits; ++i) {bit[i] = (STdata.cb.barcodes[g] >> i) & 1ULL;}
      // ... and compute gene rank as a proportion of total number of genes
      double g_prop = ((double)gi + 1.0) / N_genes;
      // For each blank ...
      for (int bi = 0; bi < N_blanks; ++bi) {
        int b = STdata.cb.blanks[bi];
        // ... compute blank rank as a proportion of total number of blanks
        double b_prop = ((double)bi + 1.0) / N_blanks;
        // Compute the relative distance between this gene and this blank's rank ... The closer they are, the more the bit-flips needed to transform g into c need to be increased.
        double rel_dist = (g_prop - b_prop) * (g_prop - b_prop);
        // Get all raw barcode reads which would be corrected to this blank
        std::vector<uint64_t> correctable_to_b = STdata.correction_table_inverted.at(b);
        // For each such raw barcode read ...
        for (uint64_t c : correctable_to_b) {
          // ... get bit-flips required to transform the true barcode into the misread barcode
          uint64_t flips_from_g_to_c = (STdata.cb.barcodes[g] ^ c) & ((1ULL << N_bits) - 1);
          int HD = __builtin_popcountll(flips_from_g_to_c);
          double w = 1.0/(rel_dist + 1.0)/(double)HD;
          std::vector<int> flip(N_bits);
          // For each bit ...
          for (int i = 0; i < N_bits; ++i) {
            // Extract whether bit i is flipped from flips_from_g_to_c
            flip[i] = (flips_from_g_to_c >> i) & 1ULL;
            if (flip[i]) {
              if (bit[i]) {
                initial_rate10[i] += w;
              } else {
                initial_rate01[i] += w;
              }
            }
            // Extract correlations between bit flips i and j, for j < i
            for (int j = 0; j < i; ++j) {
              int k = i * (i - 1) / 2 + j; // index for strict lower triangle
              if (bit[i]) {
                if (flip[i] == flip[j]) {
                  initial_corr1[k] += w;
                } else {
                  initial_corr1[k] -= w;
                }
              } else {
                if (flip[i] == flip[j]) {
                  initial_corr0[k] += w;
                } else {
                  initial_corr0[k] -= w;
                }
              }
            }
          }
        }
      }
      // Report
      if (gi % (N_genes/10) == 0 || gi == N_genes - 1) {
        Rcpp::Rcout << "  Processed gene " << gi + 1 << "/" << N_genes << "." << std::endl;
      }
    }
    // Find max values
    double max_rate10 = *std::max_element(initial_rate10.begin(), initial_rate10.end());
    double max_rate01 = *std::max_element(initial_rate01.begin(), initial_rate01.end());
    double max_corr1 = 0.0;
    double max_corr0 = 0.0;
    for (int i = 0; i < corr_free; ++i) {
      // Corr can be negative, so take absolute value when finding max
      if (std::abs(initial_corr1[i]) > max_corr1) {max_corr1 = std::abs(initial_corr1[i]);}
      if (std::abs(initial_corr0[i]) > max_corr0) {max_corr0 = std::abs(initial_corr0[i]);}
    }
    if (max_rate10 == 0.0 || max_rate01 == 0.0) {Rcpp::stop("Warning: Initial flip rate scalars are all zero. Check codebook and correction table for issues.");}
    // Normalize
    for (int i = 0; i < N_bits; ++i) {
      initial_rate10[i] = initial_rate10[i] / max_rate10;
      initial_rate01[i] = initial_rate01[i] / max_rate01;
    }
    for (int i = 0; i < corr_free; ++i) {
      if (max_corr1 > 0.0) {initial_corr1[i] = initial_corr1[i] / max_corr1;}
      if (max_corr0 > 0.0) {initial_corr0[i] = initial_corr0[i] / max_corr0;}
    }
    return {initial_rate10, initial_rate01, initial_corr1, initial_corr0};
  }

// Markov Chain Monte Carlo with Simulated Annealing to find flip rates and correlations 
FlipRates MCMCSA(
    const std::vector<double>& FR,
    const std::vector<double>& ub,
    const std::vector<double>& lb,
    const std::vector<double>& step_size,
    const std::vector<double>& temp, 
    void* data,
    int ran_seed,
    double corr_step_scale,
    double rate10_scale
  ) {
    
    // Set up steps
    int step = 0;
    int calls = 0; 
    int last_reported = -1;
    int n_steps = step_size.size(); 
    if (temp.size() != n_steps) {
      Rcpp::stop("step_size and temp vectors must be the same length");
    }
    
    // Grab data and advance sim number
    auto* d = static_cast<ST_data*>(data);
    
    // Initialize parameter vectors
    std::vector<double> FR_current = FR;
    std::vector<double> FR_next = FR;
    std::vector<double> FR_best = FR;
    int n_FR = FR.size();
    int N_bits = d->cb.N_bits;
    
    // Check that initial parameters are within bounds
    for (int i = 0; i < n_FR; ++i) {
      if (FR[i] < lb[i] || FR[i] > ub[i]) {
        Rcpp::stop("Initial parameters must be within bounds");
      }
    }
    
    // Make weight vector for msle computation, so blanks carry same overall weight as genes despite being fewer in number
    int N_barcodes = d->cb.barcodes.size();
    double blank_weight = (double)d->cb.genes.size()/(double)d->cb.blanks.size();
    std::vector<double> weights(N_barcodes, 1.0);
    for (int i : d->cb.blanks) {weights[i] = blank_weight;}
    // ... normalize 
    double weight_scaler = (double)N_barcodes/std::accumulate(weights.begin(), weights.end(), 0.0);
    for (int i = 0; i < N_barcodes; ++i) {weights[i] *= weight_scaler;}
    
    // Compute expected corrected counts from these flip rates
    FlipRates fr = pack_fr(FR_current, N_bits);
    std::vector<int> bc_counts_true = est_bc_counts_true(fr, data);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> erctc = expected_bc_counts(
      fr,
      bc_counts_true, 
      d->cb.barcodes, 
      d->correction_table_inverted,
      d->n_forks
    );
    double msle_current = compute_msle(d->bc_counts, std::get<1>(erctc), weights); 
    double msle_next = msle_current;
    double msle_least = msle_current;
    Rcpp::Rcout << "\nResampling initial parameters with MCMCSA run:\nStep: 0, msle: " << msle_current << std::endl;
    d->eval_history.msle.push_back(msle_current);
    d->eval_history.fr.push_back(FR_current);
    d->eval_history.etc.push_back(std::get<2>(erctc));
    d->eval_history.ecc.push_back(std::get<1>(erctc));
    d->eval_history.erc.push_back(std::get<0>(erctc));
    
    // Compute and save expected CR and PPV for each barcode
    d->eval_history.CR.push_back(std::vector<double>(N_barcodes, 0.0));
    d->eval_history.PPV.push_back(std::vector<double>(N_barcodes, 0.0));
    d->eval_history.est_true_bc_counts.push_back(std::vector<double>(N_barcodes, 0.0));
    for (int i = 0; i < N_barcodes; ++i) {
      d->eval_history.CR[0][i] = std::get<1>(erctc)[i] > 0.0 ? std::get<0>(erctc)[i] / std::get<1>(erctc)[i] : 0.0;
      d->eval_history.PPV[0][i] = std::get<1>(erctc)[i] > 0.0 ? std::get<2>(erctc)[i] / std::get<1>(erctc)[i] : 0.0;
      d->eval_history.est_true_bc_counts[0][i] = bc_counts_true[i];
    }
    
    // Start random-number generator and initialize a uniform distribution
    std::mt19937 rng(ran_seed);
    std::uniform_real_distribution<> unif(0.0, 1.0);
    
    while (step < n_steps) {
      
      // Generate random step (... this is the Markov chain)
      std::normal_distribution<double> norm(0.0, step_size[step]);
      FR_next = FR_current;
      for (int i = 0; i < n_FR; ++i) {
        double sz = norm(rng);
        if (i < N_bits) {sz *= rate10_scale;}
        if (i >= 2*N_bits) {sz *= corr_step_scale;}
        FR_next[i] += sz;
        // Enforce bounds ... if out of bounds, reflect back into bounds
        if (FR_next[i] < lb[i]) {
          FR_next[i] = lb[i] + (lb[i] - FR_next[i]);
          if (FR_next[i] > ub[i]) {FR_next[i] = lb[i];}
        } else if (FR_next[i] > ub[i]) {
          FR_next[i] = ub[i] - (FR_next[i] - ub[i]);
          if (FR_next[i] < lb[i]) {FR_next[i] = ub[i];}
        }
      }
      
      // Compute expected corrected counts from these flip rates
      fr = pack_fr(FR_next, N_bits);
      bc_counts_true = est_bc_counts_true(fr, data);
      std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> erctc = expected_bc_counts(
        fr,
        bc_counts_true, 
        d->cb.barcodes, 
        d->correction_table_inverted,
        d->n_forks
      );
      std::vector<double> erc = std::get<0>(erctc);
      msle_next = compute_msle(d->bc_counts, std::get<1>(erctc), weights); 
      
      // Calculate acceptance probability
      // ... idea: When msle decreases, probability of acceptance is 1; this formula
      //      controls how quickly the probability of acceptance decreases as the mse increases
      double acceptance_prob = std::min(1.0, std::exp(-(msle_next - msle_current)/temp[step]));
     
      // Accept or reject the proposed step
      if (unif(rng) < acceptance_prob) {
        // Accept the new parameters ... this is updating for the Markov chain
        FR_current = FR_next;
        // Update msle
        msle_current = msle_next;
        if (msle_current < msle_least) {
          msle_least = msle_current;
          FR_best = FR_current;
        }
        // Advance step 
        step++;
        // Save results ... this is step 2 of the Monte Carlo method: aggregate results
        d->eval_history.msle.push_back(msle_current);
        d->eval_history.fr.push_back(FR_current);
        d->eval_history.etc.push_back(std::get<2>(erctc));
        d->eval_history.ecc.push_back(std::get<1>(erctc));
        d->eval_history.erc.push_back(std::get<0>(erctc));
        // ... compute and save expected CR and PPV for each barcode
        int N_barcodes = d->cb.barcodes.size();
        d->eval_history.CR.push_back(std::vector<double>(N_barcodes, 0.0));
        d->eval_history.PPV.push_back(std::vector<double>(N_barcodes, 0.0));
        d->eval_history.est_true_bc_counts.push_back(std::vector<double>(N_barcodes, 0.0));
        for (int i = 0; i < N_barcodes; ++i) {
          d->eval_history.CR[step][i] = std::get<1>(erctc)[i] > 0.0 ? std::get<0>(erctc)[i] / std::get<1>(erctc)[i] : 0.0;
          d->eval_history.PPV[step][i] = std::get<1>(erctc)[i] > 0.0 ? std::get<2>(erctc)[i] / std::get<1>(erctc)[i] : 0.0;
          d->eval_history.est_true_bc_counts[step][i] = bc_counts_true[i];
        }
      }
      if (last_reported < step && (step % d->report_freq == 0 || step == 10)) {
        Rcpp::Rcout << "  Step: " << step << "/" << n_steps << ", msle: " << msle_current << std::endl;
        last_reported = step;
      }
      calls++;
      
    }
    Rcpp::Rcout << "\nAcceptance rate (aim for >0.2 and <0.3): " << (double)n_steps / (double)calls << std::endl;
    Rcpp::Rcout << "Best msle: " << msle_least << std::endl;
    
    // Pack best flip-rates and return as FlipRates struct
    FlipRates fr_best = pack_fr(FR_best, N_bits);
    return fr_best;
    
  }

// [[Rcpp::export]]
List mQC( 
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance,
    NumericVector step_size, 
    NumericVector temp, 
    double max_fr,
    double max_corr,
    double initial_corr,
    double corr_step_scale,
    double rate10_scale,
    int n_steps,
    int n_forks,
    int ran_seed
  ) {
    
    // Load in data
    auto STdata = load_STdata(bc_counts, codebook, max_correctable_Hamming_distance);
    // ... and hyperparameters into STdata struct for use in objective function
    STdata.n_forks = n_forks;
    STdata.report_freq = n_steps/10;
    
    // Initialize rate10 and rate01 with some reasonable starting values, e.g., 0.01 for all bits
    int N_bits = STdata.cb.N_bits; 
    int corr_free = N_bits * (N_bits - 1) / 2; // Number of relevant correlation parameters (strict lower triangle)
    auto initial_rates = scale_initial_fr(STdata); 
    size_t n = 2*N_bits + 2*corr_free; // Two sets of flip rates, two sets of bit-flip correlations
    std::vector<double> rates_plus_corr(n, 0.0); // First N_bits are rate10, second N_bits are rate01, then corr_free bits for corr1 and corr_free bits for corr0
    for (int i = 0; i < N_bits; ++i) {
      rates_plus_corr[i] = max_fr * 0.5 * std::get<0>(initial_rates)[i] * rate10_scale; // rate10, approximately 1/5th of rate01
      rates_plus_corr[N_bits + i] = max_fr * 0.5 * std::get<1>(initial_rates)[i]; // rate01
    }
    for (int i = 0; i < corr_free; ++i) {
      rates_plus_corr[2*N_bits + i] = initial_corr * std::get<2>(initial_rates)[i]; 
      rates_plus_corr[2*N_bits + corr_free + i] = initial_corr * std::get<3>(initial_rates)[i];
    }
    
    // Set bounds for parameters
    std::vector<double> lb(n, 0.0);
    std::vector<double> ub(n, 0.0); 
    for (int i = 0; i < N_bits; ++i) {
      ub[i] = max_fr * rate10_scale; 
      ub[N_bits + i] = max_fr;
    }
    for (int i = 2*N_bits; i < n; ++i) {
      lb[i] = std::max(-max_corr, -1.0); 
      ub[i] = std::min(max_corr, 1.0);  
    }
    
    // Initialize step_size and temp vectors 
    std::vector<double> step_size_schedule(n_steps, step_size[0]); 
    std::vector<double> temp_schedule(n_steps, temp[0]);
    for (int i = 1; i < n_steps; ++i) {
      step_size_schedule[i] = step_size_schedule[i - 1] + step_size[1];;
      temp_schedule[i] = temp_schedule[i - 1] + temp[1];
      if (step_size_schedule[i] < step_size[2]) {step_size_schedule[i] = step_size[2];}
      if (temp_schedule[i] < temp[2]) {temp_schedule[i] = temp[2];}
    }
    
    // Run Markov Chain Monte Carlo with Simulated Annealing to find flip rates that minimize msle between observed and expected corrected counts
    auto fr = MCMCSA(
      rates_plus_corr, 
      ub, 
      lb, 
      step_size_schedule, 
      temp_schedule, 
      &STdata,
      ran_seed,
      corr_step_scale,
      rate10_scale
    );
    
    // Check eval history length
    if (STdata.eval_history.msle.size() != n_steps + 1 ||
        STdata.eval_history.fr.size() != n_steps + 1 ||
        STdata.eval_history.ecc.size() != n_steps + 1 ||
        STdata.eval_history.erc.size() != n_steps + 1 ||
        STdata.eval_history.CR.size() != n_steps + 1 ||
        STdata.eval_history.PPV.size() != n_steps + 1 ||
        STdata.eval_history.est_true_bc_counts.size() != n_steps + 1) {
      Rcpp::stop(
        std::string("Length of eval history does not match number of steps.") +
        "\nmsle history length: " + std::to_string(STdata.eval_history.msle.size()) +
        "\nfr history length: " + std::to_string(STdata.eval_history.fr.size()) +
        "\necc history length: " + std::to_string(STdata.eval_history.ecc.size()) +
        "\nerc history length: " + std::to_string(STdata.eval_history.erc.size()) +
        "\nCR history length: " + std::to_string(STdata.eval_history.CR.size()) +
        "\nPPV history length: " + std::to_string(STdata.eval_history.PPV.size()) +
        "\nest_true_bc_counts history length: " + std::to_string(STdata.eval_history.est_true_bc_counts.size())
      ); 
    }
    
    // Copy raw data into R data structures 
    Rcpp::Rcout << "\nCollecting parameter samples." << std::endl;
    int N_barcodes = STdata.cb.barcodes.size();
    NumericMatrix FR(n_steps, n);
    NumericMatrix ecc(n_steps, N_barcodes);
    NumericMatrix erc(n_steps, N_barcodes);
    NumericMatrix CR(n_steps, N_barcodes);
    NumericMatrix PPV(n_steps, N_barcodes);
    NumericMatrix est_true_bc_counts(n_steps, N_barcodes);
    for (int i = 0; i < n_steps; ++i) {
      std::copy(
        STdata.eval_history.fr[i].begin(),
        STdata.eval_history.fr[i].end(),
        FR.row(i).begin()
      );
      std::copy(
        STdata.eval_history.ecc[i].begin(),
        STdata.eval_history.ecc[i].end(),
        ecc.row(i).begin()
      );
      std::copy(
        STdata.eval_history.erc[i].begin(),
        STdata.eval_history.erc[i].end(),
        erc.row(i).begin()
      );
      std::copy(
        STdata.eval_history.CR[i].begin(),
        STdata.eval_history.CR[i].end(),
        CR.row(i).begin()
      );
      std::copy(
        STdata.eval_history.PPV[i].begin(),
        STdata.eval_history.PPV[i].end(),
        PPV.row(i).begin()
      );
      std::copy(
        STdata.eval_history.est_true_bc_counts[i].begin(),
        STdata.eval_history.est_true_bc_counts[i].end(),
        est_true_bc_counts.row(i).begin()
      );
    }
    
    // Make data frame to hold ST data
    DataFrame STdataR = DataFrame::create(
      _["barcode"] = STdata.cb.barcodes,
      _["species"] = STdata.cb.species,
      _["barcode_rate"] = STdata.bc_rates,
      _["barcode_variance"] = STdata.bc_variance,
      _["count_observed"] = STdata.bc_counts
    );
    
    return List::create(
      _["STdata"] = STdataR,
      _["fliprates"] = FR,
      _["erc"] = erc,
      _["ecc"] = ecc,
      _["etc"] = est_true_bc_counts,
      _["CR"] = CR,
      _["PPV"] = PPV,
      _["msle"] = STdata.eval_history.msle
    );
  }

// [[Rcpp::export]]
List mQC_init(
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance,
    double max_fr,
    double initial_corr,
    double rate10_scale,
    int n_forks
  ) {
    // Heap-allocate STdata so it survives past this call
    auto* d = new ST_data(load_STdata(bc_counts, codebook, max_correctable_Hamming_distance));
    d->n_forks = n_forks;
    d->report_freq = 1;

    int N_bits    = d->cb.N_bits;
    int corr_free = N_bits * (N_bits - 1) / 2;
    auto ir       = scale_initial_fr(*d);          // (rate10, rate01, corr1, corr0)
    size_t n      = 2*N_bits + 2*corr_free;
    std::vector<double> params(n, 0.0);

    for (int i = 0; i < N_bits; ++i) {
      params[i]          = max_fr * 0.5 * std::get<0>(ir)[i] * rate10_scale;
      params[N_bits + i] = max_fr * 0.5 * std::get<1>(ir)[i];
    }
    for (int i = 0; i < corr_free; ++i) {
      params[2*N_bits + i]             = initial_corr * std::get<2>(ir)[i];
      params[2*N_bits + corr_free + i] = initial_corr * std::get<3>(ir)[i];
    }

    // XPtr takes ownership; R's GC will call `delete` when the pointer is collected
    Rcpp::XPtr<ST_data> xptr(d, true);

    return List::create(
      _["params"]   = params,
      _["data_ptr"] = xptr
    );
  }

// [[Rcpp::export]]
double mQC_msle(
    NumericVector params,
    SEXP data_ptr_sexp
  ) {
    Rcpp::XPtr<ST_data> xptr(data_ptr_sexp);
    ST_data* d     = xptr.get();
    int N_bits     = d->cb.N_bits;
    int N_barcodes = d->cb.barcodes.size();

    std::vector<double> FR(params.begin(), params.end());
    FlipRates fr = pack_fr(FR, N_bits);

    std::vector<int> bc_counts_true = est_bc_counts_true(fr, static_cast<void*>(d));
    auto erctc = expected_bc_counts(
      fr,
      bc_counts_true,
      d->cb.barcodes,
      d->correction_table_inverted,
      d->n_forks
    );

    // Uniform weights (ignoring blank upweighting)
    std::vector<double> weights(N_barcodes, 1.0);
    return compute_msle(d->bc_counts, std::get<1>(erctc), weights);
  }
