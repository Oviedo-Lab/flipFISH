
// Rcpp
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <cmath>
#include <unordered_set>
#include <unistd.h>     // fork, pipe, read, write, close
#include <sys/wait.h>   // waitpid
#include <nlopt.hpp>    // L-BFGS and other gradient-based optimizers
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

struct EvalHist {
  std::vector<double>              msle;
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
  std::vector<double> bc_rates;       // Expected count, per cell
  std::vector<double> bc_variance;    // Expected variance in count, among cells
  std::vector<int>    bc_counts;      // Vector of length N_barcodes, giving total counts for each barcode across all cells
  // Codebook
  Codebook                                       cb;
  int                                            max_correctable_Hamming_distance;
  std::unordered_map<uint64_t, int>              correction_table;
  std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted;
  // Parameter estimation
  EvalHist         eval_history;
  std::vector<int> bc_counts_true;
  double           best_msle;
  int              n_forks;
  int              report_freq;
};

struct FlipRates {
  std::vector<double>              rate10;       // P(1 -> 0) for each bit, so length = N_bits
  std::vector<double>              rate01;       // P(0 -> 1) for each bit, so length = N_bits
  std::vector<double>              log_rate10;   // Pre-computed log(rate10): log P(1->0) for each bit, used in TR
  std::vector<double>              log_rate01;   // Pre-computed log(rate01): log P(0->1) for each bit, used in TR
  std::vector<double>              corr1;        // Strict lower triangle of luminance noise correlations when bit i is 1 
  std::vector<double>              corr0;        // ... when bit i is 0
  std::vector<std::vector<double>> log_inv_corr; // Pre-computed vector of summed j<i values log(1 - corr(i,j)), used in TR calculation
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

// Parse count data and codebook into ST_data structure
ST_data load_STdata(
    NumericMatrix bc_count_data,
    IntegerMatrix codebook,
    int           max_correctable_Hamming_distance
  ) {
    
    // Load codebook as packed integers
    Codebook cb = pack_codebook(codebook);
    
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
      eval_history, std::vector<int>(N_barcodes, 0), // placeholders
      std::numeric_limits<double>::infinity(), 0, 0  // placeholders
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
    int N_bits    = fr.rate10.size();
    int corr_free = N_bits * (N_bits - 1) / 2;
    double log_tr = 0.0;
    
    // Per-bit quantities: compute log inverse correlation, stored for backward pass.
    std::vector<double> log_inv_corr(N_bits, 0.0);
    std::vector<int>    bit_vec(N_bits), flip_vec(N_bits);
    std::vector<double> log_adj_flip_vec(N_bits);
    std::vector<bool>   clamped(N_bits, false);
    
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
    const uint64_t               BOI,                       // Barcode of interest (BOI) for which we want to compute expected count after correction
    const FlipRates&             fr,                        // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>&      bc_counts_true,            // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,             // Vector giving all possible true spot barcodes
    const std::vector<uint64_t>& corrected_to_BOI,          // Vector of barcodes that would be corrected to barcode of interest (BOI)
    std::vector<double>*         grad_ecc_b       = nullptr
  ) {
    int    N_bits          = fr.rate10.size();
    int    N_barcodes      = true_barcodes.size();
    int    N_correctable   = corrected_to_BOI.size();
    double count_corrected = 0.0;
    double count_read      = 0.0;
    double count_true      = 0.0;
    // Index i ranges over barcodes corrected to BOI, index j ranges over all codebook barcodes 
    for (int j = 0; j < N_barcodes; ++j) {
      if (bc_counts_true[j] > 0) {
        double scale = (double)bc_counts_true[j];
        double tr    = 0.0;
        for (int i = 0; i < N_correctable; ++i) {
          // Get bit-flips required to transform the true barcode into the misread barcode
          uint64_t flips_to_BOI = (true_barcodes[j] ^ corrected_to_BOI[i]) & ((1ULL << N_bits) - 1);
          // TR call: also accumulates scale * ∂TR/∂params into grad_ecc_b when non-null.
          // ... grad_ecc_b tracks ∂count_corrected/∂params because count_corrected = Σ_j Σ_i scale * TR.
          double tr_ = TR(true_barcodes[j], flips_to_BOI, fr, grad_ecc_b, scale);
          tr += tr_;
          if (corrected_to_BOI[i] == BOI) {count_read += tr_ * scale;}
        }
        count_corrected += tr * scale;
        if (true_barcodes[j] == BOI) {count_true += tr * scale;}
      }
    }
    return {count_read, count_corrected, count_true};
  }

// Estimate expected barcode counts after correction, as a function of flip rates and true barcode counts, for all barcodes.
// grad_ecc (optional): pre-zeroed flat matrix stored row-major, size N_barcodes * n_params.
//   When non-null, ∂ecc[b]/∂params is accumulated into the slice (*grad_ecc)[b*n_params .. (b+1)*n_params-1].
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> expected_bc_counts(
    const FlipRates&                                      fr,                                   // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>&                               bc_counts_true,                       // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>&                          true_barcodes,                        // Vector giving all possible true spot barcodes
    const std::unordered_map<int, std::vector<uint64_t>>& correction_table_inverted,            // Inverted correction table mapping each barcode index to vector of misread barcodes that would be corrected to it
    int                                                   n_forks,
    std::vector<double>*                                  grad_ecc                   = nullptr  // if non-null: flat [N_barcodes × n_params], pre-zeroed
  ) {
    int N_barcodes = true_barcodes.size();
    int N_bits     = fr.rate10.size();
    int corr_free  = N_bits * (N_bits - 1) / 2;
    int n_params   = 2*N_bits + 2*corr_free;
    bool do_grad   = (grad_ecc != nullptr);
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
      // Run in parallel with forking.
      // Child payload layout (per barcode in batch):
      //   [erc_0..erc_{N-1} | ecc_0..ecc_{N-1} | etc_0..etc_{N-1} | grad_b0 .. grad_b{N-1}]
      // where each grad_b is n_params doubles. Gradient section present only when do_grad=true.
      
      // Pipes for inter-process communication
      std::vector<int>                pids(n_forks);
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
          int payload_size = (do_grad ? 3 + n_params : 3) * N_barcodes_batch;
          std::vector<double> erctc_child(payload_size, 0.0);
          std::vector<double> grad_b(do_grad ? n_params : 0, 0.0); // per-barcode gradient buffer
          
          for (int b = 0; b < N_barcodes_batch; ++b) {
            if (do_grad) std::fill(grad_b.begin(), grad_b.end(), 0.0);
            std::tuple<double, double, double> erctc = expected_bc_count(
              true_barcodes[barcode_batches[i][b]],
              fr, bc_counts_true, true_barcodes, 
              correction_table_inverted.at(barcode_batches[i][b]),
              do_grad ? &grad_b : nullptr
            );
            erctc_child[b]                      = std::get<0>(erctc); // expected read count
            erctc_child[b + N_barcodes_batch]   = std::get<1>(erctc); // expected corrected count
            erctc_child[b + 2*N_barcodes_batch] = std::get<2>(erctc); // expected true count
            if (do_grad) {
              std::copy(grad_b.begin(), grad_b.end(),
                        erctc_child.begin() + 3*N_barcodes_batch + b*n_params);
            }
          } 
          
          // Send result 
          const char* buffer        = reinterpret_cast<const char*>(erctc_child.data());
          size_t      nbytes        = sizeof(double) * erctc_child.size();
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
        std::vector<double> erctc_child(payload_size, 0.0);
        
        // Read the row from the pipe into the buffer
        char*  buffer     = reinterpret_cast<char*>(erctc_child.data());
        size_t nbytes     = sizeof(double) * erctc_child.size();
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
          erc[bc_idx] += erctc_child[b];
          ecc[bc_idx] += erctc_child[b + N_barcodes_batch];
          etc[bc_idx] += erctc_child[b + 2*N_barcodes_batch];
          if (do_grad) {
            std::copy(
              erctc_child.begin() + 3*N_barcodes_batch + b*n_params,
              erctc_child.begin() + 3*N_barcodes_batch + (b+1)*n_params,
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
        std::tuple<double, double, double> erctc = expected_bc_count(
          true_barcodes[b],
          fr, bc_counts_true, true_barcodes, 
          correction_table_inverted.at(b),
          do_grad ? &grad_b : nullptr
        );
        erc[b] = std::get<0>(erctc); // expected read count
        ecc[b] = std::get<1>(erctc); // expected corrected count
        etc[b] = std::get<2>(erctc); // expected true counts
        if (do_grad) {
          std::copy(grad_b.begin(), grad_b.end(), (*grad_ecc).begin() + b * n_params);
        }
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
    void*            data
  ) {
    // Grab data
    auto*       d             = static_cast<ST_data*>(data);
    const auto& true_barcodes = d->cb.barcodes;
    int         N_barcodes    = true_barcodes.size();
    int         N_bits        = d->cb.N_bits;
    
    // Find expected decoding rate
    double mean_Hamming_weight = 0.0; 
    for (int i = 0; i < N_barcodes; ++i) {mean_Hamming_weight += (double)__builtin_popcountll(true_barcodes[i]);}
    mean_Hamming_weight /= (double)N_barcodes;
    // ... find expected flip rate
    double mean_rate10 = std::accumulate(fr.rate10.begin(), fr.rate10.end(), 0.0);
    mean_rate10       /= (double)N_bits;
    double mean_rate01 = std::accumulate(fr.rate01.begin(), fr.rate01.end(), 0.0);
    mean_rate01       /= (double)N_bits;
    double expected_flip_rate = (mean_rate10*mean_Hamming_weight + mean_rate01*((double)N_bits - mean_Hamming_weight)) / (double)N_bits;
    // ... find degree of freedom
    double mean_corr1  = std::accumulate(fr.corr1.begin(), fr.corr1.end(), 0.0);
    mean_corr1        /= (double)fr.corr1.size();
    double mean_corr0  = std::accumulate(fr.corr0.begin(), fr.corr0.end(), 0.0);
    mean_corr0        /= (double)fr.corr0.size();
    double DoF         = (1.0 - mean_corr1)*mean_Hamming_weight + (1.0 - mean_corr0)*((double)N_bits - mean_Hamming_weight);
    // ... find expected decoding rate
    double expected_decoding_rate = R::ppois(
      d->max_correctable_Hamming_distance, // Max number of bits that can flip and still be decoded
      expected_flip_rate * DoF,            // Expected number of bit flips per spot, given the expected flip rate and number of bits
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

// nlopt objective callback for L-BFGS minimization of mQC_msle.
//    'f_data' must point to an ST_data struct.
//    nlopt signals "gradient not needed" by passing an empty grad vector.
static double nlopt_mQC_msle(
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
    
    FlipRates            fr              = pack_fr(x, N_bits);
    //std::vector<int>     bc_counts_true  = est_bc_counts_true(fr, static_cast<void*>(d));
    std::vector<double>  grad_ecc_flat;
    if (do_grad) grad_ecc_flat.assign(N_barcodes * n_params, 0.0);
    
    auto erctc = expected_bc_counts(
      fr,
      d->bc_counts_true,
      d->cb.barcodes,
      d->correction_table_inverted,
      d->n_forks,
      do_grad ? &grad_ecc_flat : nullptr
    );
    
    double msle = compute_msle(d->bc_counts, std::get<1>(erctc));
    
    if (do_grad) {
      std::fill(grad.begin(), grad.end(), 0.0);
      const auto& ecc = std::get<1>(erctc);
      for (int b = 0; b < N_barcodes; ++b) {
        // Chain rule: ∂MSLE/∂θ_k = (1/N_B) * 2*(log(ecc_b+1)-log(obs_b+1))/(ecc_b+1) * ∂ecc_b/∂θ_k
        double coeff = 2.0 * (std::log(ecc[b] + 1.0) - std::log((double)d->bc_counts[b] + 1.0))
                       / (ecc[b] + 1.0) / (double)N_barcodes;
        for (int k = 0; k < n_params; ++k) {
          grad[k] += coeff * grad_ecc_flat[b * n_params + k];
        }
      }
    }
    
    d->eval_history.msle.push_back(msle);
    if (msle < d->best_msle) {
      int call_n = d->eval_history.msle.size();
      if (call_n == 1 || call_n % d->report_freq == 0) {
        Rcpp::Rcout << "  eval " << call_n << ", msle: " << msle << std::endl;
      }
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
    double        max_fr,
    double        max_corr,
    int           n_forks,
    int           maxeval   = 1000,
    double        ftol_rel  = 1e-8,
    double        xtol_rel  = 1e-6
  ) {
    
    // Load in data
    auto STdata = load_STdata(bc_counts, codebook, max_correctable_Hamming_distance);
    STdata.n_forks     = n_forks;
    STdata.report_freq = 1; //std::max(1, maxeval / 10); // print ~10 progress lines
    
    // Initialize parameters: rate10 = 0.01, rate01 = 0.05, corr = 0
    int    N_bits     = STdata.cb.N_bits;
    int    corr_free  = N_bits * (N_bits - 1) / 2;
    size_t n          = 2*N_bits + 2*corr_free;
    std::vector<double> x0(n, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      x0[i]          = 0.01;
      x0[N_bits + i] = 0.05;
    }
    // corr parameters stay at 0
    
    // Estimate true counts based on initial flip rates 
    STdata.bc_counts_true = est_bc_counts_true(pack_fr(x0, N_bits), &STdata);
    
    // Parameter bounds
    std::vector<double> lb(n, 0.0);
    std::vector<double> ub(n, 0.0); 
    for (int i = 0; i < N_bits; ++i) {
      ub[i]          = max_fr;
      ub[N_bits + i] = max_fr;
    }
    for (size_t i = 2*N_bits; i < n; ++i) {
      lb[i] = std::max(-max_corr, -1.0 + std::numeric_limits<double>::epsilon());
      ub[i] = std::min( max_corr,  1.0 - std::numeric_limits<double>::epsilon());
    }
    
    // Run L-BFGS via nlopt
    Rcpp::Rcout << "\nRunning L-BFGS (nlopt::LD_LBFGS), maxeval=" << maxeval
                << ", ftol_rel=" << ftol_rel << ", xtol_rel=" << xtol_rel << std::endl;
    nlopt::opt opt(nlopt::LD_LBFGS, n);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(nlopt_mQC_msle, &STdata);
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
    
    // Compute final quantities from optimal parameters
    FlipRates fr_best = pack_fr(x0, N_bits);
    auto erctc = expected_bc_counts(
      fr_best,
      STdata.bc_counts_true,
      STdata.cb.barcodes,
      STdata.correction_table_inverted,
      STdata.n_forks
    );
    const auto& erc_final = std::get<0>(erctc);
    const auto& ecc_final = std::get<1>(erctc);
    const auto& etc_final = std::get<2>(erctc);
    
    int N_barcodes = STdata.cb.barcodes.size();
    std::vector<double> CR_final(N_barcodes), PPV_final(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      CR_final[i]  = ecc_final[i] > 0.0 ? erc_final[i] / ecc_final[i] : 0.0;
      PPV_final[i] = ecc_final[i] > 0.0 ? etc_final[i] / ecc_final[i] : 0.0;
    }
    
    // Pack results into 1-row matrices (one row = the single optimal solution)
    NumericMatrix FR_mat(1, n),  ecc_mat(1, N_barcodes), erc_mat(1, N_barcodes),
                  etc_mat(1, N_barcodes), CR_mat(1, N_barcodes), PPV_mat(1, N_barcodes);
    for (int k = 0; k < (int)n; ++k)    FR_mat(0, k)  = x0[k];
    for (int i = 0; i < N_barcodes; ++i) {
      ecc_mat(0, i) = ecc_final[i];
      erc_mat(0, i) = erc_final[i];
      etc_mat(0, i) = etc_final[i];
      CR_mat(0, i)  = CR_final[i];
      PPV_mat(0, i) = PPV_final[i];
    }
    
    DataFrame STdataR = DataFrame::create(
      _["barcode"]          = STdata.cb.barcodes,
      _["species"]          = STdata.cb.species,
      _["barcode_rate"]     = STdata.bc_rates,
      _["barcode_variance"] = STdata.bc_variance,
      _["count_observed"]   = STdata.bc_counts
    );
    
    return List::create(
      _["STdata"]    = STdataR,
      _["fliprates"] = FR_mat,
      _["erc"]       = erc_mat,
      _["ecc"]       = ecc_mat,
      _["etc"]       = etc_mat,
      _["CR"]        = CR_mat,
      _["PPV"]       = PPV_mat,
      _["msle"]      = minf
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
