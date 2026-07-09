
// Rcpp
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <cmath>
#include <unordered_set>
#include <unistd.h>     // fork, pipe, read, write, close
#include <sys/wait.h>   // waitpid
#include <nlopt.hpp>    // L-BFGS and other gradient-based optimizers
#include <boost/math/distributions/normal.hpp>
using namespace Rcpp;
using namespace Eigen;

// Hamming distance
// int dist = __builtin_popcountll(a ^ b);

/*
 * *********************************************************************************************************************
 * Basic data structures
 */

struct Codebook {
  int                      N_bits;
  std::vector<uint64_t>    barcodes;
  std::vector<std::string> species;
  std::vector<int>         blanks;
  std::vector<int>         genes; 
};

struct EvalResults {
  std::vector<double> msle_hist;
  std::vector<double> ehc;                 // expected hit count
  std::vector<double> ecc;                 // expected corrected count
  std::vector<double> erc;                 // expected read count
  std::vector<double> CR;
  std::vector<double> PPV;
  int                 n_evals = 0; 
};

struct ST_data {
  // ST data
  std::vector<double> bc_rates;            // Expected count, per cell
  std::vector<double> bc_variance;         // Expected variance in count, among cells
  std::vector<int>    bc_counts;           // Vector of length N_barcodes, giving total counts for each barcode across all cells
  // Codebook
  Codebook                                       cb;
  int                                            max_correctable_Hamming_distance;
  std::unordered_map<uint64_t, int>              correction_table;
  std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted;
  // Parameter estimation
  EvalResults         eval_results;
  std::vector<int>    bc_counts_true;      // For storing 
  double              best_msle;
  int                 n_forks;             // For computing expected counts in parallel
  int                 max_flips;           // For approximating expected counts estimate by ignoring highly unlikely misreads
  int                 report_freq;
};

struct FlipRates {
  std::vector<double>              rate10;       // P(1 -> 0) for each bit, so length = N_bits
  std::vector<double>              rate01;       // P(0 -> 1) for each bit, so length = N_bits
  std::vector<double>              log_rate10;   // Pre-computed log(rate10): log P(1->0) for each bit, used in TR
  std::vector<double>              log_rate01;   // Pre-computed log(rate01): log P(0->1) for each bit, used in TR
  std::vector<double>              corr1;        // Strict lower triangle of luminance noise correlations when bit i is 1 
  std::vector<double>              corr0;        // ... when bit i is 0
}; 

struct FlipRatePriors {
  double expected10 = 0.01;
  double expected01 = 0.05;
  double sd10 = 0.1;
  double sd01 = 0.2;
  double expectedcorr = 0.0;
  double sdcorr = 0.2;
};

struct SpotSim {
  // For each barcode, the total number of spots read as that barcode, corrected to that barcode, and correctly read as that barcode
  std::vector<int> read_counts;
  std::vector<int> corrected_counts;
  std::vector<int> hit_counts;
};

/*
 * *********************************************************************************************************************
 * Helper functions
 */

double compute_msle(
    const std::vector<int>&    obs_counts,   // observed corrected counts
    const std::vector<double>& pred_counts   // predicted corrected counts
  ) {
    int N_barcodes = obs_counts.size();
    if (pred_counts.size() != N_barcodes) {Rcpp::stop("obs_counts and pred_counts must be the same length.");}
    double msle = 0.0;
    for (int b = 0; b < N_barcodes; ++b) {
      double le = std::log(pred_counts[b] + 1.0) - std::log(static_cast<double>(obs_counts[b]) + 1.0);
      msle += le * le;
    }
    msle /= static_cast<double>(N_barcodes);
    return msle;
  }

std::vector<int> grep_idx(
    Rcpp::CharacterVector x,
    std::string           s,
    bool                  neg = false
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

FlipRates pack_fr(
    const std::vector<double>& rates_plus_corr,
    int                        N_bits
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
    // Precompute log flip rates used in TR
    fr.log_rate10 = std::vector<double>(N_bits, 0.0);
    fr.log_rate01 = std::vector<double>(N_bits, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      fr.log_rate10[i] = std::log(fr.rate10[i]);
      fr.log_rate01[i] = std::log(fr.rate01[i]);
    }
    return fr;
  }

FlipRatePriors pack_fr_priors(
    const List& fliprate_priors
  ) {
    FlipRatePriors fr_priors; 
    if (fliprate_priors.containsElementNamed("expected10"))
      fr_priors.expected10 = as<double>(fliprate_priors["expected10"]);
    if (fliprate_priors.containsElementNamed("expected01"))
      fr_priors.expected01 = as<double>(fliprate_priors["expected01"]);
    if (fliprate_priors.containsElementNamed("sd10"))
      fr_priors.sd10 = as<double>(fliprate_priors["sd10"]);
    if (fliprate_priors.containsElementNamed("sd01"))
      fr_priors.sd01 = as<double>(fliprate_priors["sd01"]);
    if (fliprate_priors.containsElementNamed("expectedcorr"))
      fr_priors.expectedcorr = as<double>(fliprate_priors["expectedcorr"]);
    if (fliprate_priors.containsElementNamed("sdcorr"))
      fr_priors.sdcorr = as<double>(fliprate_priors["sdcorr"]);
    return fr_priors; 
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
    uint64_t               x,
    int                    n_bits,
    int                    dist,
    int                    start_bit,
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
    int      n_bits,
    int      hamming_dist
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

// ... written by ChatGPT
// [[Rcpp::export]]
Eigen::MatrixXd correlation_from_params(
    const std::vector<double>& theta,
    int n
  ) {
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
    int k = 0;
    for (int i = 0; i < n; ++i){
      L(i, i) = 1.0;
      for (int j = 0; j < i; ++j) {L(i, j) = theta[k++];}
    }
    Eigen::MatrixXd S = L * L.transpose();
    Eigen::VectorXd sd = S.diagonal().cwiseSqrt();
    Eigen::MatrixXd R = S.array().colwise() / sd.array();
    R = R.array().rowwise() / sd.transpose().array();
    return R;
  }

/*
 * *********************************************************************************************************************
 * Functions for spot decoding
 */

// Extract unique Hamming distances between barcodes
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

// Convert codebook matrix into CodeBook object
Codebook pack_codebook(
    const IntegerMatrix& codebook,
    bool                 verbose = true
  ) {
    // Extract basic info
    int             N_bits     = codebook.ncol();
    int             N_barcodes = codebook.nrow();
    CharacterVector species    = rownames(codebook);
    // Initialize new Codebook
    Codebook cb;
    cb.N_bits = N_bits;
    // Make index for blanks and genes
    cb.blanks = grep_idx(species, "Blank");
    cb.genes  = grep_idx(species, "Blank", true);
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

// Quick extract unique Hamming distances from codebook
// [[Rcpp::export]]
IntegerVector unique_Hamming_cb(
    const IntegerMatrix& codebook
  ) {
    Codebook         cb                = pack_codebook(codebook, false);
    std::vector<int> hamming_distances = unique_Hamming(cb.barcodes);
    return wrap(hamming_distances);
  }

// Build table used to correct read spot barcodes
std::unordered_map<uint64_t, int> build_correction_table(
    const Codebook& cb,
    int             max_correctable_Hamming_distance
  ) {
    std::unordered_map<uint64_t, int> correction_table;
    std::unordered_set<uint64_t> ambiguous;
    for (size_t i = 0; i < cb.barcodes.size(); ++i) {
      correction_table[cb.barcodes[i]] = i; // Exact match
      // Generate all barcodes within max_correctable_Hamming_distance
      for (int d = 1; d <= max_correctable_Hamming_distance; ++d) {
        std::vector<uint64_t> n = neighbors(cb.barcodes[i], cb.N_bits, d);
        for (int j = 0; j < n.size(); ++j) {
          if (ambiguous.count(n[j])) {
            continue; // Already removed; don't re-add
          }
          if (correction_table.count(n[j])) {
            // If this neighbor is already mapped to a different barcode, we have a tie
            if (correction_table[n[j]] != (int)i) {
              correction_table.erase(n[j]);
              ambiguous.insert(n[j]);
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
 * Functions to analytically compute expected count read from flip rates and true counts
 */

// Select random flip rates and bit-flip correlations 
std::vector<double> draw_random_fr(
    const ST_data&        STdata,
    const FlipRatePriors& fr_prior
  ) {
    // Initialize normal distributions
    std::mt19937                     rng(12345);
    std::normal_distribution<double> norm10(std::sqrt(fr_prior.expected10), fr_prior.sd10);
    std::normal_distribution<double> norm01(std::sqrt(fr_prior.expected01), fr_prior.sd01);
    std::normal_distribution<double> normCorr(fr_prior.expectedcorr, fr_prior.sdcorr);
    
    // Randomly select flip rates
    int    N_bits     = STdata.cb.N_bits;
    int    corr_free  = N_bits * (N_bits - 1) / 2;
    size_t n          = 2*N_bits + 2*corr_free;
    double near_one   = 1.0 - std::numeric_limits<double>::epsilon();
    double FRsqrt;
    std::vector<double> FR(n, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      FRsqrt         = norm10(rng);
      FR[i]          = std::min(FRsqrt*FRsqrt, near_one);
      FRsqrt         = norm01(rng);
      FR[N_bits + i] = std::min(FRsqrt*FRsqrt, near_one);
    }
    // ... and random bit-flip correlations 
    std::vector<double> corr1(corr_free); 
    std::vector<double> corr0(corr_free); 
    for (int i = 0; i < corr_free; ++i) {
      corr1[i] = std::max(std::min(normCorr(rng), near_one), -near_one);
      corr0[i] = std::max(std::min(normCorr(rng), near_one), -near_one);
    }
    auto corr1_mat = correlation_from_params(corr1, N_bits); 
    auto corr0_mat = correlation_from_params(corr0, N_bits); 
    // ... fill into FR
    int k = N_bits * 2;
    for (int i = 1; i < N_bits; ++i) {
      for (int j = 0; j < i; ++j) {
        FR[k++] = corr1_mat(i, j);
      }
    }
    for (int i = 1; i < N_bits; ++i) {
      for (int j = 0; j < i; ++j) {
        FR[k++] = corr0_mat(i, j);
      }
    }
    return FR;
  }

// Parse count data and codebook into ST_data structure
ST_data load_STdata(
    NumericMatrix bc_count_data,
    IntegerMatrix codebook,
    int           max_flips,
    int           max_correctable_Hamming_distance,
    bool          verbose
  ) {
    
    // Load codebook as packed integers
    Codebook cb = pack_codebook(codebook, verbose);
    
    // Check max_flips value 
    if (max_flips < 1) {max_flips = cb.N_bits;}
    
    // Get info from bc_count_data 
    int             N_barcodes    = bc_count_data.nrow();
    int             rate_col      = -1;
    int             variance_col  = -1; 
    int             count_col     = -1;
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
    std::vector<double> bc_rates(N_barcodes);
    std::vector<double> bc_variance(N_barcodes);
    std::vector<int>    bc_counts(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      bc_rates[i]         = static_cast<double>(bc_count_data(i, rate_col));
      bc_variance[i]      = static_cast<double>(bc_count_data(i, variance_col));
      bc_counts[i]        = static_cast<int>(bc_count_data(i, count_col));
    }
    
    // Build correction table 
    std::unordered_map<uint64_t, int> correction_table = build_correction_table(cb, max_correctable_Hamming_distance);
    if (verbose) {
      Rcpp::Rcout << "Correction table built with max-correctable Hamming distance " << max_correctable_Hamming_distance << "." << std::endl;
      Rcpp::Rcout << "Correction table size: " << correction_table.size() << " correctable barcodes." << std::endl;
    }
    
    // Invert correction table
    std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted = invert_lookup_table(correction_table);
    if (verbose) {
      Rcpp::Rcout << "Correction table inverted." << std::endl;
    }
    
    // Return parsed ST data
    return {
      bc_rates, bc_variance, bc_counts,
      cb, 
      max_correctable_Hamming_distance,
      correction_table, correction_table_inverted,
      EvalResults(), std::vector<int>(N_barcodes, 0),              // placeholders
      std::numeric_limits<double>::infinity(), 
      1, max_flips, 10                                             // n_forks, max_flips, report_freq
    };
    
  }

// The probability that the sequence of flips transform_flips occurs for a spot with true barcode bc,
// given per-bit flip rates and inter-bit-flip correlations.
//
// Model: for each bit i (processed in order), the correlation-adjusted flip probability is
//   p_adj[i] = rate_flip[i] / prod_{j<i, flip[j]=1}(1 - corr(i,j))
// The no-flip probability is then 1 - p_adj[i], ensuring the per-bit probabilities sum to 1.
// The correction is applied only to the flip branch; the no-flip branch is derived from it.
// 
// ... rewritten by Claude Sonnet 4.6, so that adjustment by log_inv_corr applied only if the flip
// happens, and with the clamping guard. 
// grad (optional): flat parameter gradient vector [rate10 | rate01 | corr1 | corr0].
//   When non-null, grad_scale * ∂TR/∂params is accumulated additively into *grad.
//   The caller is responsible for pre-zeroing *grad before the first accumulating call.
double TR(
    const uint64_t       bc,
    const uint64_t       transform_flips,
    const FlipRates&     fr,
    std::vector<double>* grad             = nullptr,
    double               grad_scale       = 1.0
  ) {
    int    N_bits    = fr.rate10.size();
    int    corr_free = N_bits * (N_bits - 1) / 2;
    double log_tr    = 0.0;
    
    // Per-bit quantities: compute log inverse correlation, stored for backward pass.
    std::vector<double> log_inv_corr(N_bits, 0.0);
    std::vector<double> log_adj_flip_vec(N_bits);
    std::vector<bool>   clamped(N_bits, false);
    std::vector<int>    bit_vec(N_bits), flip_vec(N_bits);
    
    // Forward pass
    for (int i = 0; i < N_bits; ++i) {
      bit_vec[i]  = (bc              >> i) & 1ULL;
      flip_vec[i] = (transform_flips >> i) & 1ULL;
      for (int j = 0; j < i; ++j) {
        if (bit_vec[i] & flip_vec[j]) {
          log_inv_corr[i] += std::log(1.0 - fr.corr1[i * (i-1) / 2 + j]);
        } else if (flip_vec[j]) {
          log_inv_corr[i] += std::log(1.0 - fr.corr0[i * (i-1) / 2 + j]);
        }
      }
      // Log of the correlation-adjusted flip rate for bit i:
      //   log(rate_flip[i]) - sum_{j<i, flip[j]=1} log(1 - corr(i,j))
      double laf = (bit_vec[i] ? fr.log_rate10[i] : fr.log_rate01[i]) - log_inv_corr[i];
      // Guard: clamp adjusted flip probability to (0, 1) so both branches remain defined.
      //   If correlations push the adjusted rate above 1 the model is undefined; clamp to
      //   the largest representable value strictly below 1 so log1p(-p) stays finite.
      clamped[i]         = (laf >= 0.0);
      log_adj_flip_vec[i] = clamped[i] ? std::log(1.0 - std::numeric_limits<double>::epsilon()) : laf;
      if (flip_vec[i]) {
        log_tr += log_adj_flip_vec[i];
      } else {
        log_tr += std::log1p(-std::exp(log_adj_flip_vec[i]));
      }
    }
    double tr = std::exp(log_tr);
    
    // Backward pass (only when gradient is requested)
    if (grad) {
      for (int i = 0; i < N_bits; ++i) {
        if (clamped[i]) continue;  // δ_i = 0 in the clamped region
        double p_i     = std::exp(log_adj_flip_vec[i]);
        // δ_i = ∂log(TR)/∂log(p_i)
        double delta_i = flip_vec[i] ? 1.0 : (-p_i / (1.0 - p_i));
        double base    = grad_scale * tr * delta_i;
        // Gradient w.r.t. rate_flip[i]: ∂log(p_i)/∂rate_flip = 1/rate_flip
        if (bit_vec[i]) {
          (*grad)[i]          += base / fr.rate10[i];
        } else {
          (*grad)[N_bits + i] += base / fr.rate01[i];
        }
        // Gradient w.r.t. corr parameters: ∂log(p_i)/∂corr_relevant[k] = 1/(1-corr)
        for (int j = 0; j < i; ++j) {
          if (!flip_vec[j]) continue;  // corr[i,j] only appears when flip_j = 1
          int k = i * (i-1) / 2 + j;
          if (bit_vec[i]) {
            (*grad)[2*N_bits + k]             += base / (1.0 - fr.corr1[k]);
          } else {
            (*grad)[2*N_bits + corr_free + k] += base / (1.0 - fr.corr0[k]);
          }
        }
      }
    }
    
    return tr;
  }

// Function to compute the expected count for a given barcode of interest (BOI)
// grad_ecc_b (optional): pre-zeroed flat vector [rate10|rate01|corr1|corr0] of length n_params.
//   When non-null, ∂count_corrected/∂params is accumulated additively into *grad_ecc_b.
std::tuple<double, double, double> expected_bc_count(
    const uint64_t               BOI,                         // Barcode of interest (BOI) for which we want to compute expected count after correction
    const FlipRates&             fr,                          // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>&      bc_counts,                   // Vector of same length as barcodes, giving the ground-truth number of spots with each barcode
    const std::vector<uint64_t>& barcodes,                    // Vector giving all possible true spot barcodes (genes, not blanks)
    const std::vector<uint64_t>& corrected_to_BOI,            // Vector of barcodes that would be corrected to barcode of interest (BOI)
    int                          max_flips, 
    std::vector<double>*         grad_ecc_b        = nullptr
  ) {
    int    N_bits          = fr.rate10.size();
    int    N_barcodes      = barcodes.size();
    int    N_correctable   = corrected_to_BOI.size();
    double count_corrected = 0.0;
    double count_read      = 0.0;
    double count_hit       = 0.0;
    // Index i ranges over barcodes corrected to BOI, index j ranges over all codebook barcodes 
    for (int j = 0; j < N_barcodes; ++j) {
      if (bc_counts[j] > 0) {
        double scale = (double)bc_counts[j];
        double tr    = 0.0;
        for (int i = 0; i < N_correctable; ++i) {
          // Get bit-flips required to transform the true barcode into the misread barcode
          uint64_t flips_to_BOI = (barcodes[j] ^ corrected_to_BOI[i]) & ((1ULL << N_bits) - 1);
          // Check hamming distance 
          if (__builtin_popcountll(flips_to_BOI) > max_flips) {continue;}
          // TR call: also accumulates scale * ∂TR/∂params into grad_ecc_b when non-null.
          // ... grad_ecc_b tracks ∂count_corrected/∂params because count_corrected = Σ_j Σ_i scale * TR.
          double tr_ = TR(barcodes[j], flips_to_BOI, fr, grad_ecc_b, scale);
          tr += tr_;
          if (corrected_to_BOI[i] == BOI) {count_read += tr_ * scale;}
        }
        count_corrected += tr * scale;
        if (barcodes[j] == BOI) {count_hit += tr * scale;}
      }
    }
    return {count_read, count_corrected, count_hit};
  }

// Estimate expected barcode counts after correction, as a function of flip rates and true barcode counts, for all barcodes.
// grad_ecc (optional): pre-zeroed flat matrix stored row-major, size N_barcodes * n_params.
//   When grad_ecc is non-null, ∂ecc[b]/∂params is accumulated into the slice (*grad_ecc)[b*n_params .. (b+1)*n_params-1].
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> expected_bc_counts(
    const FlipRates&                                      fr,                         // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>&                               bc_counts,                  // Vector of same length as barcodes, giving the ground-truth number of spots with each barcode
    const std::vector<uint64_t>&                          barcodes,                   // Vector giving all possible true spot barcodes (genes, not blanks)
    const std::unordered_map<int, std::vector<uint64_t>>& correction_table_inverted,  // Inverted correction table mapping each barcode index to vector of misread barcodes that would be corrected to it
    int                                                   n_forks,
    int                                                   max_flips,
    std::vector<double>*                                  grad_ecc = nullptr          // if non-null: flat [N_barcodes × n_params], pre-zeroed
  ) {
    int N_barcodes = barcodes.size();
    int N_bits     = fr.rate10.size();
    int corr_free  = N_bits * (N_bits - 1) / 2;
    int n_params   = 2*N_bits + 2*corr_free;
    bool do_grad   = (grad_ecc != nullptr);
    std::vector<double> ecc(N_barcodes, 0.0); // expected corrected counts
    std::vector<double> erc(N_barcodes, 0.0); // expected read counts
    std::vector<double> ehc(N_barcodes, 0.0); // expected hit (i.e., correctly read) counts
    
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
      // Run in parallel with forking.
      // Child payload layout (per barcode in batch):
      //   [erc_0..erc_{N-1} | ecc_0..ecc_{N-1} | ehc_0..ehc_{N-1} | grad_b0 .. grad_b{N-1}]
      // where each grad_b is n_params doubles. Gradient section present only when do_grad=true.
      
      // Pipes for inter-process communication
      std::vector<int>                pids(n_forks);
      std::vector<std::array<int, 2>> pipes(n_forks); 
      
      // Initialize pipes 
      for (int i = 0; i < n_forks; ++i) {pipe(pipes[i].data());}
      
      // fork processes
      for (int i = 0; i < n_forks; i++) {
        pid_t pid              = fork();
        int   N_barcodes_batch = barcode_batches[i].size();
        
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
          int payload_size = (do_grad ? 3 + n_params : 3) * N_barcodes_batch;
          std::vector<double> erchc_child(payload_size, 0.0);
          std::vector<double> grad_b(do_grad ? n_params : 0, 0.0); // per-barcode gradient buffer
          
          for (int b = 0; b < N_barcodes_batch; ++b) {
            if (do_grad) std::fill(grad_b.begin(), grad_b.end(), 0.0);
            std::tuple<double, double, double> erchc = expected_bc_count(
              barcodes[barcode_batches[i][b]],
              fr, bc_counts, barcodes, 
              correction_table_inverted.at(barcode_batches[i][b]),
              max_flips,
              do_grad ? &grad_b : nullptr
            );
            erchc_child[b]                      = std::get<0>(erchc); // expected read count
            erchc_child[b + N_barcodes_batch]   = std::get<1>(erchc); // expected corrected count
            erchc_child[b + 2*N_barcodes_batch] = std::get<2>(erchc); // expected hit count
            if (do_grad) {
              std::copy(grad_b.begin(), grad_b.end(),
                        erchc_child.begin() + 3*N_barcodes_batch + b*n_params);
            }
          } 
          
          // Send result 
          const char* buffer        = reinterpret_cast<const char*>(erchc_child.data());
          size_t      nbytes        = sizeof(double) * erchc_child.size();
          size_t      total_written = 0;
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
          
        } else if (pid > 0) {    // parent process 
          pids[i] = pid;         // Grab child pid
          close(pipes[i][1]);    // Close write end
        } else { 
          Rcpp::stop("Fork failed!");
        } 
        
      } 
      
      // Fetch results from pipes
      for (int i = 0; i < n_forks; i++) {
        int N_barcodes_batch = barcode_batches[i].size();
        int payload_size     = (do_grad ? 3 + n_params : 3) * N_barcodes_batch;
        std::vector<double> erchc_child(payload_size, 0.0);
        
        // Read the row from the pipe into the buffer
        char*  buffer     = reinterpret_cast<char*>(erchc_child.data());
        size_t nbytes     = sizeof(double) * erchc_child.size();
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
          int bc_idx = barcode_batches[i][b];
          erc[bc_idx] += erchc_child[b];
          ecc[bc_idx] += erchc_child[b + N_barcodes_batch];
          ehc[bc_idx] += erchc_child[b + 2*N_barcodes_batch];
          if (do_grad) {
            std::copy(
              erchc_child.begin() + 3*N_barcodes_batch + b*n_params,
              erchc_child.begin() + 3*N_barcodes_batch + (b+1)*n_params,
              (*grad_ecc).begin() + bc_idx * n_params
            );
          }
        }
        close(pipes[i][0]);           // Close read end
        int status;
        waitpid(pids[i], &status, 0);
        if (WIFSIGNALED(status)) {Rcpp::Rcout << "Child killed by signal " << WTERMSIG(status) << std::endl;}
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {Rcpp::Rcout << "Child exited with code " << WEXITSTATUS(status) << std::endl;}
      }
      
    } else {
      // Run in serial
      std::vector<double> grad_b(do_grad ? n_params : 0, 0.0); // reused per-barcode gradient buffer
      for (int b : barcode_batches[0]) {
        if (do_grad) std::fill(grad_b.begin(), grad_b.end(), 0.0);
        std::tuple<double, double, double> erchc = expected_bc_count(
          barcodes[b],
          fr, bc_counts, barcodes, 
          correction_table_inverted.at(b),
          max_flips,
          do_grad ? &grad_b : nullptr
        );
        erc[b] = std::get<0>(erchc); // expected read count
        ecc[b] = std::get<1>(erchc); // expected corrected count
        ehc[b] = std::get<2>(erchc); // expected hit counts
        if (do_grad) {
          std::copy(grad_b.begin(), grad_b.end(), (*grad_ecc).begin() + b * n_params);
        }
      }
    }
    
    return {erc, ecc, ehc};
  }

// Estimate misread loss scalar 
double estimate_misread_loss(
    const ST_data& STdata,
    int            max_flips,
    const List&    fliprate_priors
  ) {
    // Grab and set counts ... uses observed corrected counts as ground-truth counts
    std::vector<int> BCcounts = STdata.bc_counts;
    for (int i : STdata.cb.blanks) {BCcounts[i] = 0;}
    int total_spots_gt = std::accumulate(BCcounts.begin(), BCcounts.end(), 0);
    // Draw and pack random flip rates 
    int       N_bits     = STdata.cb.N_bits;
    int       N_barcodes = STdata.cb.barcodes.size(); 
    auto      FR         = draw_random_fr(STdata, pack_fr_priors(fliprate_priors));
    FlipRates fr         = pack_fr(FR, N_bits);
    // Compute expected counts
    auto expected_counts = expected_bc_counts(
      fr, 
      BCcounts, 
      STdata.cb.barcodes, 
      STdata.correction_table_inverted, 
      STdata.n_forks, 
      max_flips
    );
    int total_spots_ecc = 0;
    for (int i = 0; i < N_barcodes; ++i) {
      total_spots_ecc += std::get<1>(expected_counts)[i];
    }
    return static_cast<double>(total_spots_gt) / static_cast<double>(total_spots_ecc);
  }

/*
 * *********************************************************************************************************************
 * Functions to estimate flip rates, PPV, and CR
 */

// nlopt objective callback for L-BFGS minimization of msle.
//    Note: 'f_data' must point to an ST_data struct.
//    nlopt signals "gradient not needed" by passing an empty grad vector.
static double mQC_msle(
    const std::vector<double>& x,
    std::vector<double>&       grad,
    void*                      f_data
  ) {
    ST_data* d          = static_cast<ST_data*>(f_data);
    int      N_bits     = d->cb.N_bits;
    int      N_barcodes = d->cb.barcodes.size();
    int      corr_free  = N_bits * (N_bits - 1) / 2;
    int      n_params   = 2*N_bits + 2*corr_free;
    bool     do_grad    = !grad.empty();
    
    FlipRates fr = pack_fr(x, N_bits);
    std::vector<double>  grad_ecc_flat;
    if (do_grad) grad_ecc_flat.assign(N_barcodes * n_params, 0.0);
    
    auto erchc = expected_bc_counts(
      fr,
      d->bc_counts_true,
      d->cb.barcodes,
      d->correction_table_inverted,
      d->n_forks,
      d->max_flips,
      do_grad ? &grad_ecc_flat : nullptr
    );
    
    double msle = compute_msle(d->bc_counts, std::get<1>(erchc));
    
    if (do_grad) {
      std::fill(grad.begin(), grad.end(), 0.0);
      const auto& ecc = std::get<1>(erchc);
      for (int b = 0; b < N_barcodes; ++b) {
        // Chain rule: ∂MSLE/∂θ_k = (1/N_B) * 2*(log(ecc_b+1)-log(obs_b+1))/(ecc_b+1) * ∂ecc_b/∂θ_k
        double coeff = 2.0 * (std::log(ecc[b] + 1.0) - std::log((double)d->bc_counts[b] + 1.0))
                       / (ecc[b] + 1.0) / (double)N_barcodes;
        for (int k = 0; k < n_params; ++k) {
          grad[k] += coeff * grad_ecc_flat[b * n_params + k];
        }
      }
    }
    
    // Advance eval counter
    d->eval_results.n_evals++; 
    if (msle < d->best_msle) {
      int call_n = d->eval_results.n_evals;
      if (call_n == 1 || call_n % d->report_freq == 0) {
        Rcpp::Rcout << "  eval " << call_n << ", msle: " << msle << std::endl;
      }
      // Save eval history 
      d->eval_results.msle_hist.push_back(msle);
      d->eval_results.ehc = std::get<2>(erchc);
      d->eval_results.ecc = std::get<1>(erchc);
      d->eval_results.erc = std::get<0>(erchc);
      // Compute and save expected CR and PPV for each barcode
      d->eval_results.CR  = std::vector<double>(N_barcodes, 0.0);
      d->eval_results.PPV = std::vector<double>(N_barcodes, 0.0);
      for (int i = 0; i < N_barcodes; ++i) {
        d->eval_results.CR[i]  = std::get<1>(erchc)[i] > 0.0 ? std::get<0>(erchc)[i] / std::get<1>(erchc)[i] : 0.0;
        d->eval_results.PPV[i] = std::get<1>(erchc)[i] > 0.0 ? std::get<2>(erchc)[i] / std::get<1>(erchc)[i] : 0.0;
      }
      // Save best_msle
      d->best_msle = msle;
    }
    
    return msle;
  }

// ... revised by Claude Sonnet 4.6
// [[Rcpp::export]]
List mQC( 
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int           max_correctable_Hamming_distance,
    int           n_forks,
    int           max_flips, 
    int           report_freq = 1,
    int           maxeval     = 1000,
    List          fliprate_priors = List()
  ) {
    
    // Load in data
    auto STdata = load_STdata(
      bc_counts, 
      codebook, 
      max_flips,
      max_correctable_Hamming_distance,
      true // Make verbose printouts
      );
    STdata.n_forks     = n_forks;
    STdata.report_freq = report_freq;
    
    // Initialize parameters: rate10 = 0.01, rate01 = 0.05, corr = 0
    int N_bits     = STdata.cb.N_bits;
    int N_barcodes = STdata.cb.barcodes.size(); 
    int corr_free  = N_bits * (N_bits - 1) / 2;
    size_t n       = 2*N_bits + 2*corr_free;
    std::vector<double> x0(n, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      x0[i]          = 0.01;
      x0[N_bits + i] = 0.05;
    }
    std::vector<double> x0_ = x0; 
    
    // Parameter bounds
    double near_one = 1.0 - std::numeric_limits<double>::epsilon(); 
    std::vector<double> lb(n, 0.0);
    std::vector<double> ub(n, 0.0); 
    for (int i = 0; i < N_bits; ++i) {
      ub[i]          =  near_one;
      ub[N_bits + i] =  near_one;
    }
    for (size_t i = 2*N_bits; i < n; ++i) {
      lb[i]          = -near_one;
      ub[i]          =  near_one;
    }
    
    // Estimate misread loss scaler 
    double misread_loss_scaler = estimate_misread_loss(STdata, max_flips, fliprate_priors);
    
    // Rescale observed counts by flip rates to compensate for misread loss
    for (int i = 0; i < N_barcodes; ++i) {
      STdata.bc_counts_true[i] = static_cast<int>(std::round(static_cast<double>(STdata.bc_counts[i]) * misread_loss_scaler));
    }
    // ... set blanks to zero
    for (int idx : STdata.cb.blanks) {
      STdata.bc_counts_true[idx] = 0;
    }
    
    // Initialize matrices to hold results 
    NumericMatrix erchc_plus(N_barcodes, 5);
    colnames(erchc_plus) = CharacterVector({"erc", "ecc", "ehc", "CR", "PPV"});
    NumericVector fr(n);
    
    // Run L-BFGS via nlopt
    double ftol_rel = 1e-8;
    double xtol_rel = 1e-6;
    Rcpp::Rcout << "\nRunning L-BFGS (nlopt::LD_LBFGS)" << ", maxeval=" << maxeval 
                << ", ftol_rel=" << ftol_rel << ", xtol_rel=" << xtol_rel << std::endl;
    nlopt::opt opt(nlopt::LD_LBFGS, n);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(mQC_msle, &STdata);
    opt.set_ftol_rel(ftol_rel);
    opt.set_xtol_rel(xtol_rel);
    opt.set_maxeval(maxeval);
    
    double minf = std::numeric_limits<double>::quiet_NaN();
    try {
      nlopt::result res = opt.optimize(x0, minf);
      Rcpp::Rcout << "\nL-BFGS finished (result code " << (int)res << "), final msle: " << minf << std::endl;
    } catch (const std::exception& e) {
      Rcpp::Rcout << "\nL-BFGS warning: " << e.what()
                  << "\nProceeding with best parameters found so far." << std::endl;
    }
   
    // Pack results
    for (int k = 0; k < (int)n; ++k) {fr(k)  = x0[k];}
    for (int i = 0; i < N_barcodes; ++i) {
      erchc_plus(i, 0) = STdata.eval_results.erc[i];
      erchc_plus(i, 1) = STdata.eval_results.ecc[i];
      erchc_plus(i, 2) = STdata.eval_results.ehc[i];
      erchc_plus(i, 3) = STdata.eval_results.CR[i];
      erchc_plus(i, 4) = STdata.eval_results.PPV[i];
    }
    
    return List::create(
      _["STdata"]             = DataFrame::create(
        _["barcode"]          = STdata.cb.barcodes,
        _["species"]          = STdata.cb.species,
        _["barcode_rate"]     = STdata.bc_rates,
        _["barcode_variance"] = STdata.bc_variance,
        _["count_observed"]   = STdata.bc_counts
      ),
      _["fliprates"]          = fr,
      _["erchc_plus"]         = erchc_plus,
      _["msle"]               = minf
    );
  }

/*
 * *********************************************************************************************************************
 * Correlation matrix functions
 */

// Multivariate normal CDF, upper tail
double mvnorm_cdf_uppertail(
    const NumericVector& threshold, 
    const NumericMatrix& sigma    // covariance matrix of dimension n less than 1000
  ) {
    
    if (sigma.nrow() != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (sigma.nrow() >= 1000) {Rcpp::stop("Covariance matrix must be less than 1000x1000");}
    if (sigma.nrow() != threshold.size()) {Rcpp::stop("Matrix diagonal and threshold vector must be same length");}
    
    Function pmvnorm("pmvnorm", Environment::namespace_env("mvtnorm"));
    // ... uses Genz algorithm
    
    double prob = as<double>(
      pmvnorm( // by default, lower = -Inf, upper = Inf, and mean = 0.
        Named("lower") = threshold, 
        Named("sigma") = sigma, 
        Named("keepAttr") = false
      )
    );
    
    return prob;
    
  }

// Normal CDF, with inverse
double norm_cdf(
    double x,
    double mu,
    double sd,
    bool   inverse
  ) {
    double xc = x;
    using boost::math::normal; 
    normal standard_normal(mu, sd);
    if (inverse) {
      if (xc < 1e-10) {xc = 1e-10;}
      if (xc > 1 - 1e-10) {xc = 1.0 - 1e-10;}
      return boost::math::quantile(standard_normal, xc);
    } else {
      return boost::math::cdf(standard_normal, xc);
    }
  }

// For estimating sigma for dichotomized Gaussian simulation
NumericVector dg_sigma_formula(
    double               threshold, // threshold for dichotomization
    const NumericVector& cov,       // desired covarance after dichotomization
    const NumericMatrix& sigma      // covariance matrix of multivariate Gaussian
  ) {
    // We know threshold and cov. By finding the sigma which sends this function 
    //   to zero, we can find the covariance needed for dichotomized Gaussian simulation
   
    // Check dimension
    int dim  = sigma.nrow();
    if (dim != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (dim != cov.size())   {Rcpp::stop("Covariance vector must have the same length as sigma diagonal");}
    
    // Find probability of a point being above the threshold along all dimensions
    double Phi2_upper = mvnorm_cdf_uppertail(
      Rcpp::rep(threshold, dim), 
      sigma
    );
   
    // Find probability of a point being below the threshold along one dimension
    double Phi = norm_cdf(
      threshold, 
      0.0,     // mean
      1.0,     // sd
      false    // return inverse? No, return cdf
    );
    
    // Desired sigma will be the one which sends all elements to zero
    //  Formula: cov = Phi2_upper - (1 - Phi) * (1 - Phi) is derived as follows: 
    //   By definition of cov, cov = E[X1*X2] - E[X1]*E[X2].
    //   In this case, E[X1] = E[X2] = P(X > threshold) = 1 - Phi.
    //   X1*X2 != 0 only if both X1 and X2 > threshold, which occurs with probability Phi2_upper.
    NumericVector residuals(dim);
    for (int i = 0; i < dim; i++) {
      residuals[i] = cov[i] - Phi2_upper + (1.0 - Phi) * (1.0 - Phi);
    }
    
    return residuals;
    
  }

// Wrapper for use with find-root-bisection algorithm 
double dg_sigma_formula_scalar(
    double threshold,               // threshold for dichotomization
    double cov,                     // desired covarance after dichotomization
    double sigma                    // Gaussian covariance
  ) {
    
    // Construct covariance matrix sigma (2x2)
    NumericMatrix sigmaMat(2, 2);
    sigmaMat(_,0) = NumericVector::create(1.0, sigma);
    sigmaMat(_,1) = NumericVector::create(sigma, 1.0);
    
    // Include self-covariance on front, which is the variance, sd^2
    //   The Gaussian is normal, so mean is zero and sd is 1, so variance is 1
    NumericVector cov1 = {1.0, cov};
    
    // Evaluate formula and return second value
    NumericVector residual = dg_sigma_formula(threshold, cov1, sigmaMat);
    // Return only the second element, corresponding to cov
    return residual[1]; 
    
  }

// Function to find sigma by root bisection 
double dg_find_sigma_RootBisection(
    double threshold,               // threshold for dichotomization
    double cov                      // desired covarance after dichotomization
  ) {
    
    // Set search parameters 
    const int max_iter = 50; 
    const double tol = 1e-4;
    
    // Initiate sigmas
    double sigma_lower = -0.999;
    double sigma_upper = 0.999;
    
    // Evaluate formula
    double fx_lower = dg_sigma_formula_scalar(threshold, cov, sigma_lower);
    double fx_upper = dg_sigma_formula_scalar(threshold, cov, sigma_upper);
    
    // Run checks 
    if (abs(fx_lower) < tol) {return sigma_lower;}
    else if (abs(fx_upper) < tol) {return sigma_upper;}
    else if (fx_lower * fx_upper > tol) {return 0.0;} // Both initial covariance values lie on same side of zero crossing
    
    // Run bisection
    double fx = std::numeric_limits<double>::infinity();
    double sigma_mid;
    int iter = 0;
    while (abs(fx) > tol && iter < max_iter) {
      
      // Find midpoint
      sigma_mid = (sigma_lower + sigma_upper)/2.0;
      fx = dg_sigma_formula_scalar(threshold, cov, sigma_mid);
      
      // Update bounds
      if (fx > 0.0) {
        sigma_lower = sigma_mid;
      } else {
        sigma_upper = sigma_mid;
      }
      
      // Update iteration
      iter++;
      
    }
    
    return sigma_mid; 
    
  }


/*
 * *********************************************************************************************************************
 * DG simulation functions and related stuff
 */

std::vector<int> simulate_spots_for_barcode_b(
    int                                      b,                     // ID of barcode to simulate
    int                                      count,                 // Number of spots to simulate
    const MatrixXd&                          noised_corr1_Cholesky,
    const MatrixXd&                          noised_corr0_Cholesky,
    const std::unordered_map<uint64_t, int>& correction_table, 
    const std::vector<uint64_t>&             barcodes, 
    std::mt19937&                            rng
  ) {
    
    // Barcode info
    const int N_barcodes = barcodes.size();
    const int N_bits     = noised_corr1_Cholesky.cols();
    
    // Vector to hold read, corrected, and hit counts for each barcode
    std::vector<int> rch_counts(3 * N_barcodes, 0); 
    const int        read_offset      = 0;
    const int        corrected_offset = N_barcodes;
    const int        hit_offset       = 2 * N_barcodes;
    
    if (count >= 1) {
      
      // Sample luminance levels from multivariate normal 
      // ... take standard normal samples
      int d = N_bits;
      std::vector<double> Z_flat(count*d);
      std::normal_distribution<double> norm(0.0, 1.0);
      for (int i = 0; i < count*d; ++i) {Z_flat[i] = norm(rng);}
      // ... correlate + shift mean
      std::vector<std::vector<double>> lum(N_bits, std::vector<double>(count)); // outer vector (cols) as bits, index of inner (rows) as spots
      for (int j = 0; j < N_bits; ++j) {
        // Extract bit j of barcode b
        int bit = (barcodes[b] >> j) & 1ULL;
        // ... for each spot i with true identity of barcode b ...
        for (int i = 0; i < count; ++i) {
          // Set luminance of bit j for spot i according to expected bit value: 0 -> -1, 1 ->  1
          lum[j][i] += (double)(bit * 2.0 - 1.0);
          // Apply noise in accordance with the expected flip rate for bit j and bit-flip correlations
          for (int k = 0; k < N_bits; ++k) {
            int bit_k = (barcodes[b] >> k) & 1ULL;
            if (bit_k) {
              lum[j][i] += Z_flat[i * N_bits + k] * noised_corr1_Cholesky(k, j);
            } else {
              lum[j][i] += Z_flat[i * N_bits + k] * noised_corr0_Cholesky(k, j);
            }
          }
          
        }
      }
      
      // Simulate the spots for this barcode
      for (int k = 0; k < count; ++k) {
        
        // Decode luminance values and pack into barcode integer
        uint64_t spot_bc = 0;
        for (int j = 0; j < N_bits; ++j) {spot_bc |= (uint64_t)(lum[j][k] > 0.0) << j;}
        uint64_t spot_bc_corrected = spot_bc;
        int label_corrected;
        int label_read = -1;
        
        // Correct decoded label
        auto it = correction_table.find(spot_bc);
        if (it == correction_table.end()) {label_corrected = -1;} // uncorrectable
        else {label_corrected = it->second;}
        
        // ... and correct decoded barcode
        if (label_corrected >= 0) {
          spot_bc_corrected = barcodes[label_corrected];
          if (spot_bc_corrected == spot_bc) {label_read = label_corrected;}
        }
        
        // Add (accumulate) labels
        if (label_read >= 0) {++rch_counts[read_offset + label_read];}
        if (label_corrected >= 0) {
          ++rch_counts[corrected_offset + label_corrected];
          if (label_corrected == b) {++rch_counts[hit_offset + label_corrected];}
        }
        
      }
    }
    
    return rch_counts;
  }

SpotSim make_SpotSim(
    const std::vector<int>&                  bc_counts,             // Ground-truth counts, per barcode
    const Codebook&                          cb, 
    const FlipRates&                         fr, 
    const std::unordered_map<uint64_t, int>& correction_table,
    int                                      ran_seed,
    int                                      n_forks
  ) {
    
    // Initialize spot sim
    SpotSim sim;
    
    // Extract hyperparameters 
    int total_spots = std::accumulate(bc_counts.begin(), bc_counts.end(), 0);
    int N_barcodes  = cb.barcodes.size();
    int N_bits      = cb.N_bits;
    
    // Make bit_noise matrix 
    // ... compute inverse of target flip rates using normal-distribution quantiles (inverse CDF)
    Rcpp::Rcout << "Computing inverse of target flip rates using normal-distribution quantiles (inverse CDF)." << std::endl;
    std::vector<double> rate10_inv(N_bits);
    std::vector<double> rate01_inv(N_bits);
    for (int i = 0; i < N_bits; ++i) {
      rate10_inv[i] = R::qnorm(fr.rate10[i], 0.0, 1.0, 1, 0); // qnorm(p, mean = 0, sd = 1, lower_tail = true, log_p = false)
      rate01_inv[i] = R::qnorm(1 - fr.rate01[i], 0.0, 1.0, 1, 0); 
    }
    // ... set bit noise by barcode
    Rcpp::Rcout << "Setting bit noise by barcode." << std::endl; 
    std::vector<std::vector<double>> bit_noise(N_barcodes, std::vector<double>(N_bits));
    for (int i = 0; i < N_barcodes; ++i) {
      for (int j = 0; j < N_bits; ++j) {
        // Extract bit ... barcode i, bit j
        int bit = (cb.barcodes[i] >> j) & 1ULL;
        // Convert: 0 -> -1, 1 ->  1
        double m = (double)bit * 2.0 - 1.0;
        if (bit == 1) {
          bit_noise[i][j] = -m / rate10_inv[j];
        } else {  
          bit_noise[i][j] = -m / rate01_inv[j];
        }
      }
    }
    
    // Construct batches of barcodes to simulate in parallel
    std::vector<std::vector<int>> barcode_batches(n_forks);
    for (int i = 0; i < N_barcodes; ++i) {barcode_batches[i % n_forks].push_back(i);}
    if (n_forks == 1 && barcode_batches[0].size() != N_barcodes) {
      Rcpp::stop("Error in batching barcodes for simulation. Expected all barcodes to be in one batch when n_forks=1.");
    }
    
    // Build noised covariance matrices and their Cholesky decompositions for each barcode
    // STOPPED HERE! NOT CORRECT. Can't just plug fr.corr1 in. At least, must be explicit: the correlation between luminance noise is not the same as the correlation in bit-flips!! As written above, fr.corr1 and fr.corr0 are not luminance noise correlation (must correct that comment), they are bit-flip correlations. Must invert the relationship I proved earlier to get the continuous-valued luminance noise correlations which would produce the stipulated bit-flip correlations. 
    Rcpp::Rcout << "Building noised covariance matrices and their Cholesky decomposition for each barcode." << std::endl; 
    std::vector<MatrixXd> noised_corr1_Cholesky_per_barcode(N_barcodes);
    std::vector<MatrixXd> noised_corr0_Cholesky_per_barcode(N_barcodes);
    MatrixXd              noised_corr1(N_bits, N_bits);
    MatrixXd              noised_corr0(N_bits, N_bits);
    for (int b = 0; b < N_barcodes; ++b) {
      for (int i = 0; i < N_bits; ++i) {
        for (int j = 0; j <= i; ++j) {
          if (i == j) {
            noised_corr1(i, j) = bit_noise[b][i] * bit_noise[b][i];
            noised_corr0(i, j) = bit_noise[b][i] * bit_noise[b][i];
          } else {
            int k = i * (i - 1) / 2 + j;
            noised_corr1(i, j) = bit_noise[b][i] * fr.corr1[k] * bit_noise[b][j];
            noised_corr1(j, i) = noised_corr1(i, j);
            noised_corr0(i, j) = bit_noise[b][i] * fr.corr0[k] * bit_noise[b][j];
            noised_corr0(j, i) = noised_corr0(i, j);
          }
        }
      }
      // Cholesky decomposition
      Eigen::LLT<MatrixXd> llt1(noised_corr1);
      Eigen::LLT<MatrixXd> llt0(noised_corr0);
      // ... check positive-definiteness
      if (llt1.info() != Eigen::Success) {
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(noised_corr1);
        double min_eval = es.eigenvalues().minCoeff();
        Rcpp::Rcout << "Cholesky failure\n";
        Rcpp::Rcout << "barcode: " << b << "\n";
        Rcpp::Rcout << "min eigenvalue: " << min_eval << "\n" << std::endl;
        Rcpp::stop(
          "Cholesky failed in mvn sampling."
        );
      }
      if (llt0.info() != Eigen::Success) {
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(noised_corr0);
        double min_eval = es.eigenvalues().minCoeff();
        Rcpp::Rcout << "Cholesky failure\n";
        Rcpp::Rcout << "barcode: " << b << "\n";
        Rcpp::Rcout << "min eigenvalue: " << min_eval << "\n" << std::endl;
        Rcpp::stop(
          "Cholesky failed in mvn sampling."
        );
      }
      noised_corr1_Cholesky_per_barcode[b] = llt1.matrixL();
      noised_corr0_Cholesky_per_barcode[b] = llt0.matrixL();
    }
    
    // For each barcode, simulate and decode spot counts
    Rcpp::Rcout << "Simulating and decoding spot counts for each barcode." << std::endl; 
    int cache_size = 3 * N_barcodes;
    std::vector<int> rch_counts(cache_size, 0);
    if (n_forks > 1) {
      // Run in parallel with forking
      
      // Pipes for inter-process communication
      std::vector<int>                pids(n_forks);
      std::vector<std::array<int, 2>> pipes(n_forks); 
      
      // Initialize pipes 
      for (int i = 0; i < n_forks; ++i) {pipe(pipes[i].data());}
      
      // fork processes
      for (int i = 0; i < n_forks; i++) {
        pid_t pid = fork();
        
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
          std::vector<int> rch_counts_child(cache_size, 0);
          
          for (int b : barcode_batches[i]) {
            std::mt19937 rng(ran_seed + i*n_forks + b);
            std::vector<int> temp_vec = simulate_spots_for_barcode_b(
              b, 
              bc_counts[b], 
              noised_corr1_Cholesky_per_barcode[b],
              noised_corr0_Cholesky_per_barcode[b],
              correction_table, 
              cb.barcodes, 
              rng
            ); 
            for (int j = 0; j < cache_size; ++j) {rch_counts_child[j] += temp_vec[j];}
          }
          
          // Send result 
          const char* buffer = reinterpret_cast<const char*>(rch_counts_child.data());
          size_t nbytes = sizeof(int) * cache_size;
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
        std::vector<int> temp_vec(cache_size, 0);
        
        // Read the row from the pipe into the buffer
        char* buffer = reinterpret_cast<char*>(temp_vec.data());
        size_t nbytes = sizeof(int) * cache_size;
        size_t total_read = 0;
        while (total_read < nbytes) {
          ssize_t n_read = read(
            pipes[i][0],
            buffer + total_read,
            nbytes - total_read
          );
          if (n_read == 0) {
            close(pipes[i][0]);
            Rcpp::stop(
              "Unexpected EOF while reading pipe. "
              "Read " + std::to_string(total_read) +
                " of " + std::to_string(nbytes) + " bytes."
            );
          }
          if (n_read <= 0) {
            close(pipes[i][0]);
            Rcpp::stop("Pipe read failed");
          }
          total_read += static_cast<size_t>(n_read);
        }
        
        for (int j = 0; j < cache_size; ++j) {rch_counts[j] += temp_vec[j];}
        close(pipes[i][0]);           // Close read end
        int status;
        waitpid(pids[i], &status, 0);
        if (WIFSIGNALED(status)) {Rcpp::Rcout << "Child killed by signal " << WTERMSIG(status) << std::endl;}
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {Rcpp::Rcout << "Child exited with code " << WEXITSTATUS(status) << std::endl;}
      }
      
    } else {
      // Run in serial
      for (int b : barcode_batches[0]) {
        std::mt19937 rng(ran_seed + b);
        std::vector<int> temp_vec = simulate_spots_for_barcode_b(
          b, 
          bc_counts[b], 
          noised_corr1_Cholesky_per_barcode[b],
          noised_corr0_Cholesky_per_barcode[b],
          correction_table, 
          cb.barcodes, 
          rng
        ); 
        for (int j = 0; j < cache_size; ++j) {rch_counts[j] += temp_vec[j];}
      }
    }
    
    // Resize vectors to count spots per barcode
    sim.read_counts.resize(N_barcodes, 0);
    sim.corrected_counts.resize(N_barcodes, 0);
    sim.hit_counts.resize(N_barcodes, 0);
    
    // Parse out accumulations to spot sim
    const int read_offset      = 0;
    const int corrected_offset = N_barcodes;
    const int hit_offset      = 2 * N_barcodes;
    for (int b = 0; b < N_barcodes; b++) {
      sim.read_counts[b]      = rch_counts[read_offset + b];
      sim.corrected_counts[b] = rch_counts[corrected_offset + b];
      sim.hit_counts[b]       = rch_counts[hit_offset + b];
    }
    
    return sim;
  }

// [[Rcpp::export]]
List test_fr_recovery(
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int           n_sims, 
    int           max_correctable_Hamming_distance,
    int           n_forks,
    int           max_flips, 
    int           report_freq,
    int           maxeval,
    List          fliprate_priors = List()
  ) {
    
    // Load in data
    auto STdata = load_STdata(
      bc_counts, 
      codebook, 
      max_flips,
      max_correctable_Hamming_distance,
      false  // Suppress verbose printouts
    );
    STdata.n_forks     = n_forks;
    STdata.report_freq = report_freq;
    
    // Find bc_counts columns 
    int             rate_col      = -1;
    int             variance_col  = -1; 
    int             count_col     = -1;
    CharacterVector data_colnames = colnames(bc_counts);
    for (int i = 0; i < bc_counts.ncol(); ++i) {
      if (     data_colnames[i] == "rates")    {rate_col     = i;} 
      else if (data_colnames[i] == "variance") {variance_col = i;} 
      else if (data_colnames[i] == "counts")   {count_col    = i;} 
    }
    // ... save original count column 
    NumericVector bc_counts_original(bc_counts(_, count_col)); 
    
    // Grab and set counts ... uses observed corrected counts as ground-truth counts
    std::vector<int> BCcounts = STdata.bc_counts;
    for (int i : STdata.cb.blanks) {BCcounts[i] = 0;}
    
    // Draw and pack random flip rates 
    int       N_bits     = STdata.cb.N_bits;
    int       N_barcodes = STdata.cb.barcodes.size(); 
    int       corr_free  = N_bits * (N_bits - 1) / 2;
    size_t    n          = 2*N_bits + 2*corr_free;
    auto      FR         = draw_random_fr(STdata, pack_fr_priors(fliprate_priors));
    FlipRates fr         = pack_fr(FR, N_bits);
    
    // Compute ground truth 
    auto ground_truth = expected_bc_counts(
      fr, 
      BCcounts, 
      STdata.cb.barcodes, 
      STdata.correction_table_inverted, 
      n_forks, 
      max_flips
    );
    std::vector<double> ecc_expected(N_barcodes, 0.0); 
    std::vector<double> PPV_expected(N_barcodes, 0.0);
    for (int i = 0; i < N_barcodes; ++i) {
      ecc_expected[i]  = std::get<1>(ground_truth)[i];
      PPV_expected[i]  = std::get<1>(ground_truth)[i] > 0.0 ? std::get<2>(ground_truth)[i] / std::get<1>(ground_truth)[i] : 0.0;
    }
    
    // Initialize matrices to hold results 
    NumericMatrix fr_est(n_sims, n);
    NumericMatrix sim_cnt(n_sims, N_barcodes);
    NumericMatrix ecc_est(n_sims, N_barcodes); 
    NumericMatrix PPV_est(n_sims, N_barcodes); 
    
    for (int s = 0; s < n_sims; ++s) {
      
      // Simulate corrected reads from this data and stipulated set of flip rates 
      int s_plus = s + 1; 
      Rcpp::Rcout << "\n------ Running Pass Number " << s_plus << "/" << n_sims << " ------" << std::endl;
      Rcpp::Rcout << "\nMaking dichotomized-Gaussian simulation with given data ..." << std::endl; 
      int ran_seed = 123;
      SpotSim sim = make_SpotSim(
        BCcounts, 
        STdata.cb, 
        fr, 
        STdata.correction_table, 
        ran_seed + s, 
        n_forks
      );
      
      // Use simulation results to rewrite bc_counts
      for (int i = 0; i < N_barcodes; ++i) {
        sim_cnt(s, i) = sim.corrected_counts[i];
        bc_counts(i, rate_col) /= bc_counts_original(i) / static_cast<double>(std::max(1, sim.corrected_counts[i]));
        bc_counts(i, count_col) = static_cast<double>(sim.corrected_counts[i]);
      }
      
      // Run mQC to recover flip rates 
      Rcpp::Rcout << "\nRunning misread QC with L-BFGS (nlopt) to recover flip rates...\n" << std::endl; 
      List res = mQC(
        bc_counts, codebook,
        max_correctable_Hamming_distance,
        n_forks, max_flips, 
        report_freq, maxeval,
        fliprate_priors
      );
      // ... extract flip-rate vector (size n)
      fr_est(s,_)              = Rcpp::as<NumericVector>(res["fliprates"]);
      // ... and ecc and PPV vectors (size N_barcodes)
      NumericMatrix erchc_plus = res["erchc_plus"];
      ecc_est(s,_)             = NumericVector(erchc_plus(_, 1));
      PPV_est(s,_)             = NumericVector(erchc_plus(_, 4));
      
      // ecc_est, PPV_est, erchc_plus are computed analytically from flip rates and ground-truth counts
      
      // Test 1: ... Want colMeans(ecc_est) ... eh, no, colMeans(sim_cnt) ... to approach ecc_expected (and same with colMeans(PPV_est) and PPV_expected)
      //   sim_cnt depends only on the flip rates and the stochasticity of the simulation, while ecc_est also depends on the model fit of flip rates
      // Test 2: ... Want colMeans(fr_est) to approach FR
      
    }
    
    return List::create(
      _["fr_est"]            = fr_est,
      _["fr_stipulated"]     = FR,
      _["PPV_est"]           = PPV_est,
      _["PPV_expected"]      = PPV_expected,  // "Expected", given the stipulated flip rates
      _["ecc_est"]           = ecc_est,       // These are expected corrected counts based on the flip rates fit to sim_cnt
      _["ecc_expected"]      = ecc_expected,  // "Expected", given the stipulated flip rates
      _["sim_counts"]        = sim_cnt        // These are (corrected) counts simulated directly from stipulated flip rates
    );
    
  }

/*
 * *********************************************************************************************************************
 * Test helpers (exported for use in testthat only)
 */

// Sum TR(bc, flips, fr) over all 2^N_bits flip patterns.
// Under the corrected model this must equal 1.0 for every barcode and every valid parameter set.
// Arguments mirror pack_fr: rate10 and rate01 are per-bit flip rates; corr1 and corr0 are the
// strict lower-triangle correlation vectors (length N_bits*(N_bits-1)/2 each).
// ... written by Claude Sonnet 4.6
// [[Rcpp::export]]
double tr_sum_check(
    int                        bc,
    const std::vector<double>& rate10,
    const std::vector<double>& rate01,
    const std::vector<double>& corr1,
    const std::vector<double>& corr0
  ) {
    int N_bits = rate10.size();
    std::vector<double> params;
    params.insert(params.end(), rate10.begin(), rate10.end());
    params.insert(params.end(), rate01.begin(), rate01.end());
    params.insert(params.end(), corr1.begin(),  corr1.end());
    params.insert(params.end(), corr0.begin(),  corr0.end());
    FlipRates fr = pack_fr(params, N_bits);
   
    double total = 0.0;
    uint64_t n_patterns = 1ULL << N_bits;
    for (uint64_t flips = 0; flips < n_patterns; ++flips) {
      total += TR(static_cast<uint64_t>(bc), flips, fr);
    }
    return total;
  }
