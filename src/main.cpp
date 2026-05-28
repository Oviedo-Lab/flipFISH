
// Rcpp
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <nlopt.hpp>
#include <cmath>
#include <unistd.h>     // fork, pipe, read, write, close
#include <sys/wait.h>   // waitpid
using namespace Rcpp;
using namespace Eigen;

// Hamming distance
// int dist = __builtin_popcountll(a ^ b);

struct Codebook {
  int N_bits;
  std::vector<uint64_t> barcodes;
  std::vector<std::string> species;
  std::vector<int> blanks;
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
  std::vector<double> eval_history;
  int n_reruns;
  int least_bc_n;
  int n_forks;
  int sim_num;
  int call_freq;
};

struct FlipRates {
  std::vector<double> rate10; // P(1 -> 0) for each bit, so length = N_bits
  std::vector<double> rate01; // P(0 -> 1) for each bit, so length = N_bits
  std::vector<double> log_rate10; // Pre-computed vector of log(rate10) values, used in logTR calculation
  std::vector<double> log_rate01; // Pre-computed vector of log(rate01) values, used in logTR calculation
  std::vector<double> log_rate11; // Pre-computed vector of 1-log(rate10) values, used in logTR calculation
  std::vector<double> log_rate00; // Pre-computed vector of 1-log(rate01) values, used in logTR calculation
  MatrixXd corr; // N_bits x N_bits matrix of luminance noise correlations between bits
  std::vector<double> log_inv_corr; // Pre-computed vector of summed j<i values log(1 - corr(i,j)), used in logTR calculation
}; 

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

FlipRates pack_fr(
    const std::vector<double>& rates_plus_corr,
    int N_bits
  ) {
    FlipRates fr;
    fr.rate10.assign(rates_plus_corr.begin(), rates_plus_corr.begin() + N_bits);
    fr.rate01.assign(rates_plus_corr.begin() + N_bits, rates_plus_corr.begin() + 2*N_bits);
    fr.corr = make_corr_matrix(
      std::vector<double>(rates_plus_corr.begin() + 2*N_bits, rates_plus_corr.end()), 
      N_bits
    );
    fr.log_inv_corr = std::vector<double>(N_bits, 0.0);
    fr.log_rate00 = std::vector<double>(N_bits, 0.0);
    fr.log_rate01 = std::vector<double>(N_bits, 0.0);
    fr.log_rate10 = std::vector<double>(N_bits, 0.0);
    fr.log_rate11 = std::vector<double>(N_bits, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      fr.log_rate10[i] = std::log(fr.rate10[i]);
      fr.log_rate01[i] = std::log(fr.rate01[i]);
      fr.log_rate00[i] = std::log(1.0 - fr.rate01[i]);
      fr.log_rate11[i] = std::log(1.0 - fr.rate10[i]);
      for (int j = 0; j < i; ++j) {
        fr.log_inv_corr[i] += std::log(1.0 - fr.corr(i,j));
      }
    }
    return fr;
  }

struct SpotSim {
  // For each barcode, the total number of spots read as that barcode, corrected to that barcode, and correctly read as that barcode
  std::vector<int> read_counts;
  std::vector<int> corrected_counts;
  std::vector<int> true_counts;
};

std::vector<int> grep_idx(
    Rcpp::CharacterVector x,
    std::string s
  ) {
    std::vector<int> idx;
    for (int i = 0; i < x.size(); ++i) {
      std::string xi = Rcpp::as<std::string>(x[i]);
      if (xi.find(s) != std::string::npos) {idx.push_back(i);}
    }
    return idx;
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

std::vector<int> simulate_spots_for_barcode_b(
    int b,                 // ID of barcode to simulate
    int count,             // Number of spots to simulate
    const MatrixXd& noised_corr_Cholesky,
    const std::unordered_map<uint64_t, int>& correction_table, 
    const std::vector<uint64_t>& barcodes, 
    std::mt19937& rng
  ) {
    
    // Barcode info
    const int N_barcodes = barcodes.size();
    const int N_bits = noised_corr_Cholesky.cols();
    
    // Vector to hold read, corrected, and true counts for each barcode
    std::vector<int> rct_counts(3 * N_barcodes, 0); 
    const int read_offset = 0;
    const int corrected_offset = N_barcodes;
    const int true_offset = 2 * N_barcodes;
    
    if (count >= 1) {
      
      // Sample luminance levels from multivariate normal 
      // ... take standard normal samples
      int n = count;
      int d = N_bits;
      std::vector<double> Z_flat(n*d);
      std::normal_distribution<double> norm(0.0, 1.0);
      for (int i = 0; i < n*d; ++i) {Z_flat[i] = norm(rng);}
      // ... correlate + shift mean
      std::vector<std::vector<double>> lum(N_bits, std::vector<double>(count)); // outer vector (cols) as bits, index of inner (rows) as spots
      for (int j = 0; j < N_bits; ++j) {
        for (int i = 0; i < count; ++i) {
          for (int k = 0; k < N_bits; ++k) {
            lum[j][i] += Z_flat[i * N_bits + k] * noised_corr_Cholesky(k, j);
          }
          // Extract bit 
          int bit = (barcodes[b] >> j) & 1ULL;
          // Convert into bit mean: 0 -> -1, 1 ->  1
          lum[j][i] += (double)(bit * 2.0 - 1.0);
        }
      }
      
      // Simulate the spots for this barcode
      for (int k = 0; k < count; ++k) {
        
        // Decode luminance values and pack into barcode integer
        uint64_t spot_bc = 0;
        for (int b = 0; b < N_bits; ++b) {spot_bc |= (uint64_t)(lum[b][k] > 0.0) << b;}
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
        if (label_read >= 0) {++rct_counts[read_offset + label_read];}
        if (label_corrected >= 0) {
          ++rct_counts[corrected_offset + label_corrected];
          if (label_corrected == b) {++rct_counts[true_offset + label_corrected];}
        }
        
      }
    }
    
    return rct_counts;
  }

SpotSim make_SpotSim(
    const std::vector<int>& bc_counts, // Ground-truth counts, per barcode
    const Codebook& cb, 
    const FlipRates& fr, 
    const std::unordered_map<uint64_t, int>& correction_table,
    int ran_seed,
    int n_forks
  ) {
    
    // Initialize spot sim
    SpotSim sim;
    
    // Extract hyperparameters 
    int total_spots = std::accumulate(bc_counts.begin(), bc_counts.end(), 0);
    int N_barcodes = cb.barcodes.size();
    int N_bits = cb.N_bits;
    
    // Make bit_noise matrix 
    // ... compute inverse of target flip rates using normal-distribution quantiles (inverse CDF)
    std::vector<double> rate10_inv(N_bits);
    std::vector<double> rate01_inv(N_bits);
    for (int i = 0; i < N_bits; ++i) {
      rate10_inv[i] = R::qnorm(fr.rate10[i], 0.0, 1.0, 1, 0); // qnorm(p, mean = 0, sd = 1, lower_tail = true, log_p = false)
      rate01_inv[i] = R::qnorm(1 - fr.rate01[i], 0.0, 1.0, 1, 0); 
    }
    // ... set bit noise by barcode
    std::vector<std::vector<double>> bit_noise(N_barcodes, std::vector<double>(N_bits));
    for (int i = 0; i < N_barcodes; ++i) {
      for (int j = 0; j < N_bits; ++j) {
        // Extract bit ... barcode i, bit b
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
    std::vector<MatrixXd> noised_corr_Cholesky_per_barcode(N_barcodes);
    MatrixXd noised_corr(N_bits, N_bits);
    for (int b = 0; b < N_barcodes; ++b) {
      for (int i = 0; i < N_bits; ++i) {
        for (int j = 0; j < N_bits; ++j) {
          noised_corr(i, j) =  bit_noise[b][i] * fr.corr(i, j) * bit_noise[b][j];
        }
      }
      // Regularize diagonal
      noised_corr.diagonal().array() += 1e-8;
      // Cholesky decomposition
      Eigen::LLT<MatrixXd> llt(noised_corr);
      // ... check positive-definiteness
      if (llt.info() != Eigen::Success) {
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(noised_corr);
        double min_eval = es.eigenvalues().minCoeff();
        Rcpp::Rcout << "Cholesky failure\n";
        Rcpp::Rcout << "barcode: " << b << "\n";
        Rcpp::Rcout << "min eigenvalue: " << min_eval << "\n" << std::endl;
        Rcpp::stop(
          "Cholesky failed in mvn sampling."
        );
      }
      noised_corr_Cholesky_per_barcode[b] = llt.matrixL();
    }
    
    // For each barcode, simulate and decode spot counts
    int cache_size = 3 * N_barcodes;
    std::vector<int> rct_counts(cache_size, 0);
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
          std::vector<int> rct_counts_child(cache_size, 0);
          
          for (int b : barcode_batches[i]) {
            std::mt19937 rng(ran_seed + i*n_forks + b);
            std::vector<int> temp_vec = simulate_spots_for_barcode_b(
              b, 
              bc_counts[b], 
              noised_corr_Cholesky_per_barcode[b],
              correction_table, 
              cb.barcodes, 
              rng
            ); 
            for (int j = 0; j < cache_size; ++j) {rct_counts_child[j] += temp_vec[j];}
          }
          
          // Send result 
          const char* buffer = reinterpret_cast<const char*>(rct_counts_child.data());
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
        
        for (int j = 0; j < cache_size; ++j) {rct_counts[j] += temp_vec[j];}
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
          noised_corr_Cholesky_per_barcode[b],
          correction_table, 
          cb.barcodes, 
          rng
        ); 
        for (int j = 0; j < cache_size; ++j) {rct_counts[j] += temp_vec[j];}
      }
    }
    
    // Resize vectors to count spots per barcode
    sim.read_counts.resize(N_barcodes, 0);
    sim.corrected_counts.resize(N_barcodes, 0);
    sim.true_counts.resize(N_barcodes, 0);
    
    // Parse out accumulations to spot sim
    const int read_offset = 0;
    const int corrected_offset = N_barcodes;
    const int true_offset = 2 * N_barcodes;
    for (int b = 0; b < N_barcodes; b++) {
      sim.read_counts[b] = rct_counts[read_offset + b];
      sim.corrected_counts[b] = rct_counts[corrected_offset + b];
      sim.true_counts[b] = rct_counts[true_offset + b];
    }
    
    return sim;
  }

std::pair<std::vector<double>, std::vector<double>> compute_CRPPV(
    const SpotSim& sim
  ) {
    // Compute CR and PPV
    int N_barcodes = sim.read_counts.size();
    std::vector<double> CR(N_barcodes, 0.0);
    std::vector<double> PPV(N_barcodes, 0.0);
    for (size_t i = 0; i < N_barcodes; ++i) {
      PPV[i] = sim.corrected_counts[i] > 0 ? double(sim.true_counts[i]) / double(sim.corrected_counts[i]) : 0.0;
      CR[i] = sim.corrected_counts[i] > 0 ? double(sim.read_counts[i]) / double(sim.corrected_counts[i]) : 0.0;
    }
    return {CR, PPV};
  }

double compute_sim_nll(
    int N_cells,
    const std::vector<double>& bc_rates,          // observed rates
    const std::vector<double>& bc_dispersion,     // observe dispersion (gamma) factors
    const std::vector<int>& sim_corrected_counts  // simulated corrected counts 
  ) {
    // Compute negative log likelihood of bc_counts, given corrected_counts as predicted counts
    double log_lik = 0.0;
    for (int i = 0; i < bc_rates.size(); ++i) {
      double pred_rate = (double)sim_corrected_counts[i] / (double)N_cells;
      if (pred_rate == 0.0) {pred_rate = 1e-12;} // Avoid zero predicted rates, which cause issues for the likelihood calculation
      double gamma_variance = pred_rate + bc_dispersion[i] * pred_rate * pred_rate;
      log_lik += std::log(poisson_gamma_integral(bc_rates[i], pred_rate, gamma_variance));
    }
    return -log_lik;
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

double compute_msle(
    const std::vector<int>& obs_counts,
    const std::vector<double>& pred_counts           // predicted corrected counts
  ) {
    int N_barcodes = obs_counts.size();
    if (pred_counts.size() != N_barcodes) {Rcpp::stop("obs_counts and pred_counts must be the same length.");}
    double msle = 0.0;
    for (int b = 0; b < N_barcodes; ++b) {
      double sle = std::log(pred_counts[b] + 1.0) - std::log((double)obs_counts[b] + 1.0);
      msle += sle * sle;
    }
    msle /= (double)N_barcodes;
    return msle;
  }

uint64_t pack(
    std::vector<int> bits
  ) {
    uint64_t packed = 0;
    for (int i = 0; i < bits.size(); ++i) {packed |= (uint64_t(bits[i]) << i);} 
    return packed;
  }

Codebook pack_codebook(
    const IntegerMatrix& codebook
  ) {
    int N_bits = codebook.ncol();
    int N_barcodes = codebook.nrow();
    CharacterVector species = rownames(codebook);
    Codebook cb;
    cb.N_bits = N_bits;
    cb.blanks = grep_idx(species, "Blank");
    cb.barcodes.reserve(N_barcodes);
    cb.species.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      std::vector<int> bits(N_bits);
      for (int j = 0; j < N_bits; ++j) {bits[j] = codebook(i, j);}
      cb.barcodes.push_back(pack(bits));
      cb.species.push_back(as<std::string>(species[i]));
    }
    return cb;
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

std::pair<std::vector<int>, std::vector<double>> make_ground_truth_counts(
    int N_cells, 
    const ST_data& STdata,
    const FlipRates& fr
  ) {
    // Extract ST data
    Codebook cb = STdata.cb;
    std::vector<double> bc_rates = STdata.bc_rates;
    std::vector<double> bc_variance = STdata.bc_variance;
    int N_barcodes = cb.barcodes.size();
    int N_bits = cb.N_bits;
    int max_correctable_Hamming_distance = STdata.max_correctable_Hamming_distance;
    // ... find expected decoding rate
    double mean_Hamming_weight = 0.0; 
    for (int i = 0; i < N_barcodes; ++i) {mean_Hamming_weight += (double)__builtin_popcountll(cb.barcodes[i]);}
    mean_Hamming_weight /= (double)N_barcodes;
    double mean_rate10 = std::accumulate(fr.rate10.begin(), fr.rate10.end(), 0.0);
    mean_rate10 /= (double)N_bits;
    double mean_rate01 = std::accumulate(fr.rate01.begin(), fr.rate01.end(), 0.0);
    mean_rate01 /= (double)N_bits;
    double expected_flip_rate = (mean_rate10*mean_Hamming_weight + mean_rate01*((double)N_bits - mean_Hamming_weight)) / (double)N_bits;
    double expected_decoding_rate = R::ppois(max_correctable_Hamming_distance, expected_flip_rate * (double)N_bits, true, false);
    // ... make bc_counts and bc_dispersion vectors
    std::vector<int> bc_counts(N_barcodes, 0);
    std::vector<double> bc_dispersion(N_barcodes, 0.0);
    double N_cells_adjusted = (double)N_cells / expected_decoding_rate;
    for (int i = 0; i < N_barcodes; ++i) {
      bc_counts[i] = (int)std::round(bc_rates[i] * N_cells_adjusted);
      bc_dispersion[i] = compute_gamma_dispersion(bc_rates[i], bc_variance[i]);
    }
    // ... Set blanks to zero
    for (int idx : cb.blanks) {
      bc_counts[idx] = 0;
    }
    return {bc_counts, bc_dispersion};
  }

double observed_counts_nll_by_sim(
    const std::vector<double>& rates_plus_corr, 
    std::vector<double>& grad,
    void* data                    
  ) {
   
    // Extract data
    auto* d = static_cast<ST_data*>(data);
    const auto& cb = d->cb;
    const auto& bc_rates = d->bc_rates;
    const auto& N_bits = cb.N_bits;
    // ... extract model hyperparameters
    int n_reruns = d->n_reruns;
    int least_bc_n = d->least_bc_n;
    int n_forks = d->n_forks;
    int sim_num = ++d->sim_num;
    int call_freq = d->call_freq;
    
    // Compute number of cells needed for a transcript count of least_bc_n for the least abundant barcode
    int N_cells = (int)std::round((double)least_bc_n / *std::min_element(bc_rates.begin(), bc_rates.end()));
   
    // Extract rate data
    FlipRates fr = pack_fr(rates_plus_corr, N_bits);
    
    // Extract simulated ground-truth counts and dispersion
    const auto& gtc = make_ground_truth_counts(N_cells, *d, fr);
    
    // Compute nll of observed data given these rates, with a simulation 
    double nll = 0.0;
    SpotSim sim;
    for (int i = 0; i < n_reruns; ++i) {
      // Run simulation with these rates
      sim = make_SpotSim(gtc.first, cb, fr, d->correction_table, 999 + i, n_forks);
      // Compute nll 
      nll += compute_sim_nll(N_cells, bc_rates, gtc.second, sim.corrected_counts);
    }
    nll /= (double)n_reruns; // Average nll across simulations to reduce noise
    
    if (sim_num % call_freq == 0) {
      Rcpp::Rcout << "Call: " << sim_num << ", nll: " << nll << std::endl;
    }
    
    // Save evaluation record
    d->eval_history.push_back(nll);
    
    // Return nll
    return nll; 
  }

using ObjectiveFn = double (*)(
  const std::vector<double>&,
  std::vector<double>&,
  void*
);

/*
 * Functions to analytically compute expected count read from flip rates and true counts
 */

// Log of the rate at which the sequence of flips transform_flips can be expected to occur for a spot with true barcode bc, given bit-flip rates rate10 and rate01 and correlation corr between bit-flips
double logTR(
    const uint64_t bc,
    const uint64_t transform_flips,
    const FlipRates& fr
  ) {
    // Get number of bits
    int N_bits = fr.rate10.size();
    // Initialize variable to hold log of transformation rate
    double log_tr = 0.0;
    for (int i = 0; i < N_bits; ++i) {
      // Extract bit i from barcode bc
      int bit = (bc >> i) & 1ULL;
      // Extract whether bit i is flipped from transform_flips
      int flip = (transform_flips >> i) & 1ULL;
      // Compute independent transformation rate for this bit 
      double log_tr_nocorr = fr.log_rate00[i];
      if (flip) {
        if (bit) {
          log_tr_nocorr = fr.log_rate10[i];
        } else {
          log_tr_nocorr = fr.log_rate01[i];
        }
      } else if (bit) {
        log_tr_nocorr = fr.log_rate11[i];
      }
      // Add adjusted log transformation rate for this bit to total log transformation rate
      log_tr += log_tr_nocorr - fr.log_inv_corr[i];
    }
    return log_tr; 
  }

// Number of B' incorrectly read as something corrected to B
double expected_bc_count(
    const FlipRates& fr,                            // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>& bc_counts_true,         // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,     // Vector giving all possible true spot barcodes
    const std::vector<uint64_t>& corrected_to_BOI   // vector of barcodes that would be corrected to barcode of interest (BOI)
  ) {
    int N_bits = fr.rate10.size();
    int N_correctable = corrected_to_BOI.size();
    int N_barcodes = true_barcodes.size();
    double count = 0.0;
    // For each possible barcode misread that would be corrected to the barcode indexed by k ...
    for (int j = 0; j < N_barcodes; ++j) {
      if (bc_counts_true[j] > 0) {
        double tr = 0.0;
        for (int i = 0; i < N_correctable; ++i) {
          // Get bit-flips required to transform the true barcode into the misread barcode
          uint64_t flips_to_BOI = (true_barcodes[j] ^ corrected_to_BOI[i]) & ((1ULL << N_bits) - 1);
          // Find expected count of spot misreads corrected to barcode k, from true barcode j, and add to total misread count
          tr += std::exp(logTR(true_barcodes[j], flips_to_BOI, fr));
        }
        count += tr * (double)bc_counts_true[j];
      }
    }
    return count;
  }

// Estimate expected barcode counts after correction,  as a function of flip rates and true barcode counts
std::vector<double> expected_bc_counts(
    const FlipRates& fr,                                                             // FlipRates struct holding rate10, rate01, and corr
    const std::vector<int>& bc_counts_true,                                          // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,                                      // Vector giving all possible true spot barcodes
    const std::unordered_map<int, std::vector<uint64_t>>& correction_table_inverted,  // Inverted correction table mapping each corrected barcode index to vector of misread barcodes that would be corrected to it
    int n_forks
  ) {
    int N_barcodes = true_barcodes.size();
    std::vector<double> expected_counts(N_barcodes, 0.0);
    int max_count = *std::max_element(bc_counts_true.begin(), bc_counts_true.end());
    
    // Construct batches of barcodes to simulate in parallel
    std::vector<std::vector<int>> barcode_batches(n_forks);
    for (int i = 0; i < N_barcodes; ++i) {barcode_batches[i % n_forks].push_back(i);}
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
          std::vector<double> expected_counts_child(N_barcodes_batch, 0.0);
          
          for (int b = 0; b < N_barcodes_batch; ++b) {
            expected_counts_child[b] = expected_bc_count(
              fr,
              bc_counts_true, 
              true_barcodes, 
              correction_table_inverted.at(barcode_batches[i][b])
            );
          } 
          
          // Send result 
          const char* buffer = reinterpret_cast<const char*>(expected_counts_child.data());
          size_t nbytes = sizeof(double) * expected_counts_child.size();
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
        std::vector<double> expected_counts_child(N_barcodes_batch, 0.0);
        
        // Read the row from the pipe into the buffer
        char* buffer = reinterpret_cast<char*>(expected_counts_child.data());
        size_t nbytes = sizeof(double) * expected_counts_child.size();
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
              "\nUnexpected EOF while reading pipe " + std::to_string(i) + " of " + std::to_string(n_forks) + ". "
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
        
        for (int b = 0; b < N_barcodes_batch; ++b) {expected_counts[barcode_batches[i][b]] += expected_counts_child[b];}
        close(pipes[i][0]);           // Close read end
        int status;
        waitpid(pids[i], &status, 0);
        if (WIFSIGNALED(status)) {Rcpp::Rcout << "Child killed by signal " << WTERMSIG(status) << std::endl;}
        if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {Rcpp::Rcout << "Child exited with code " << WEXITSTATUS(status) << std::endl;}
      }
      
    } else {
      // Run in serial
      for (int b : barcode_batches[0]) {
        expected_counts[b] = expected_bc_count(
          fr,
          bc_counts_true, 
          true_barcodes, 
          correction_table_inverted.at(b)
        );
      }
    }
    
    return expected_counts;
  }

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
    double mean_rate10 = std::accumulate(fr.rate10.begin(), fr.rate10.end(), 0.0);
    mean_rate10 /= (double)N_bits;
    double mean_rate01 = std::accumulate(fr.rate01.begin(), fr.rate01.end(), 0.0);
    mean_rate01 /= (double)N_bits;
    double expected_flip_rate = (mean_rate10*mean_Hamming_weight + mean_rate01*((double)N_bits - mean_Hamming_weight)) / (double)N_bits;
    double expected_decoding_rate = R::ppois(d->max_correctable_Hamming_distance, expected_flip_rate * (double)N_bits, true, false);
    
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

double observed_counts_nll_analytic(
    const std::vector<double>& rates_plus_corr, 
    std::vector<double>& grad,
    void* data                    
  ) {
   
    // Grab data and advance sim number
    auto* d = static_cast<ST_data*>(data);
    int sim_num = ++d->sim_num;
   
    // Extract rate data
    FlipRates fr = pack_fr(rates_plus_corr, d->cb.N_bits);
   
    // Compute expected corrected counts from these flip rates
    std::vector<double> ecc = expected_bc_counts(
      fr,
      est_bc_counts_true(fr, data), 
      d->cb.barcodes, 
      d->correction_table_inverted,
      d->n_forks
    );
    
    // Compute negative log likelihood of observed barcodes given expected corrected counts as mean of Poisson distribution
    //double nll = compute_nll(*d, ecc);
    double nll = compute_msle(d->bc_counts, ecc); 
    
    if (sim_num % d->call_freq == 0 || sim_num == 10) {
      Rcpp::Rcout << "Call: " << sim_num << ", nll: " << nll << std::endl;
    }
    
    // Save evaluation record
    d->eval_history.push_back(nll);
    
    // ... and return
    return nll;
  }

// Optimize rate10 and rate01 to minimize distance between expected_bc_counts and observed_barcodes
FlipRates estimate_flip_rates(
    ST_data& STdata,
    ObjectiveFn fn,
    double max_fr,
    double ctol,
    int max_evals,
    std::string algorithm_name
  ) {
    
    Rcpp::Rcout << "Estimating flip-rates ..." << std::endl;
   
    // Initialize rate10 and rate01 with some reasonable starting values, e.g., 0.01 for all bits
    int N_bits = STdata.cb.N_bits; 
    int corr_free = N_bits*(N_bits + 1) / 2; // Number of free parameters in the correlation matrix (symmetric plus diagonal)
    size_t n = N_bits*2 + corr_free;
    std::vector<double> rates_plus_corr(n, 0.0); // First N_bits are rate10, second N_bits are rate01, rest of bits are corr matrix
    for (int i = 0; i < N_bits; ++i) {
      rates_plus_corr[i] = 0.01; // rate10
      rates_plus_corr[N_bits + i] = 0.05; // rate01
    }
    Rcpp::Rcout << "n parameters to optimize: " << n << std::endl;
    
    // Get requested optimization algorithm
    static const std::unordered_map<std::string, nlopt::algorithm> alg_map = {
      {"COBYLA", nlopt::LN_COBYLA},
      {"CRS2", nlopt::GN_CRS2_LM},
      {"NELDERMEAD", nlopt::LN_NELDERMEAD},
      {"SBPLX", nlopt::LN_SBPLX},
      {"PRAXIS", nlopt::LN_PRAXIS},
      {"ESCH", nlopt::GN_ESCH},
      {"ISRES", nlopt::GN_ISRES},
      {"DIRECT", nlopt::GN_DIRECT},
      {"DIRECT_L", nlopt::GN_DIRECT_L}
    };
    auto oa = alg_map.find(algorithm_name);
    if (oa == alg_map.end()) {Rcpp::stop("Unknown optimization algorithm");}
    
    // Set up NLopt optimizer
    nlopt::srand(12345); // Set random seed for reproducibility of stochastic optimization algorithms
    nlopt::opt opt(oa->second, n); 
    opt.set_min_objective(fn, &STdata);
    opt.set_ftol_rel(ctol);       // stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(max_evals);   // Maximum number of evaluations to try
    if (algorithm_name == "CRS2" || algorithm_name == "ISRES") {
      Rcpp::Rcout << "Setting initial seed-population size for " << algorithm_name << " to " << n + 1 << std::endl;
      opt.set_population(n + 1);
    }
    // ... add bounds 
    std::vector<double> lb(n, 0.0);
    std::vector<double> ub(n, max_fr); 
    for (int i = 2*N_bits; i < n; ++i) {
      lb[i] = -1.0; 
      ub[i] = 1.0;  
    }
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    
    // Fit model
    int success_code = 0;
    double min_fx;
    try {
      nlopt::result sc = opt.optimize(rates_plus_corr, min_fx);
      success_code = static_cast<int>(sc);
    } catch (std::exception& e) {
      Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
      success_code = 0;
    }  
    
    // Extract rate10, rate01, and corr from rates_plus_corr vector and return
    FlipRates fr = pack_fr(rates_plus_corr, N_bits);
    
    return fr;
  }

ST_data load_STdata(
    NumericMatrix bc_count_data,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance
  ) {
   
    // Load codebook as packed integers
    Rcpp::Rcout << "\nLoading data..." << std::endl;
    Codebook cb = pack_codebook(codebook);
    Rcpp::Rcout << "Codebook loaded: " << cb.barcodes.size() << " barcodes, " << cb.N_bits << " bits, " << cb.blanks.size() << " blanks." << std::endl;
    
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
    Rcpp::Rcout << "Barcode counts loaded: " << bc_counts.size() << " barcodes." << std::endl;
    
    // Build correction table 
    std::unordered_map<uint64_t, int> correction_table = build_correction_table(cb, max_correctable_Hamming_distance);
    Rcpp::Rcout << "Correction table built with max Hamming distance " << max_correctable_Hamming_distance << "." << std::endl;
    Rcpp::Rcout << "Correction table size: " << correction_table.size() << " entries." << std::endl;
    
    // Invert correction table
    std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted = invert_lookup_table(correction_table);
    Rcpp::Rcout << "Correction table inverted." << std::endl;
   
    // Return parsed ST data
    Rcpp::Rcout << "Data loading complete.\n" << std::endl;
    std::vector<double> eval_history;
    return {
      bc_rates, bc_variance, bc_counts, 
      cb, 
      max_correctable_Hamming_distance,
      correction_table, correction_table_inverted,
      eval_history, 0, 0, 0, 0 // placeholders
      };
    
  }

// [[Rcpp::export]]
List mQC( 
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance,
    double max_fr = 0.25,
    double ctol = 1e-3,
    int max_evals = 500,
    std::string algorithm_name = "COBYLA",
    int n_reruns = 1,                           // Number of times to rerun a simulation when predicting counts as a way to estimate likelihood of observed counts given a set of bit-flip rates; values >1 help to smooth out stochastic noise
    int least_bc_n = 1,                         // Expected number of transcripts for least abundant barcode, used to set counts in simulations from rates
    int n_resamples = 100,                      // Number of times to simulate a given set of flip rates, to get distribution of CR and PPV values
    int n_forks = 4                             // Maximum number of parallel processes to fork when simulating spots; set to 1 to disable forking and run in serial
  ) {
    
    // Load in data
    auto STdata = load_STdata(bc_counts, codebook, max_correctable_Hamming_distance);
    // ... and hyperparameters into STdata struct for use in objective function
    STdata.n_reruns = n_reruns;
    STdata.least_bc_n = least_bc_n;
    STdata.n_forks = n_forks;
    STdata.sim_num = 0; 
    STdata.call_freq = max_evals/10;
   
    // Make initial estimate of flip rates, assuming no correlations or over-dispersion, to use as starting point for DG simulations
    auto fr = estimate_flip_rates(
      STdata, 
      observed_counts_nll_analytic, //observed_counts_nll_by_sim,
      max_fr,
      ctol,
      max_evals,
      algorithm_name
    );
    
    // Compute expected corrected counts from these flip rates
    std::vector<double> ecc = expected_bc_counts(
      fr,
      est_bc_counts_true(fr, &STdata), 
      STdata.cb.barcodes, 
      STdata.correction_table_inverted,
      n_forks
    );
    
    // // Compute CR and PPV values via simulation with these flip rates
    // Rcpp::Rcout << "\nSimulating with estimated flip rates to estimate CR and PPV for each barcode ..." << std::endl;
    // int N_cells = (int)std::round(10.0 / *std::min_element(STdata.bc_rates.begin(), STdata.bc_rates.end()));
    // NumericVector counts_simulated = NumericVector(STdata.bc_rates.begin(), STdata.bc_rates.end()) * (double)N_cells;
    // SpotSim sim;
    // NumericMatrix corrected_counts_sim(n_resamples, STdata.cb.barcodes.size());
    // NumericMatrix CR_sim(n_resamples, STdata.cb.barcodes.size());
    // NumericMatrix PPV_sim(n_resamples, STdata.cb.barcodes.size());
    // NumericVector nll_sim(n_resamples);
    // // Extract simulated ground-truth counts and dispersion
    // const auto& gtc = make_ground_truth_counts(N_cells, STdata, fr);
    // 
    // // Estiamte CR and PPV for these rates, across multiple simulations to get a distribution of values
    // std::pair<std::vector<double>, std::vector<double>> CRPPV;
    // for (int i = 0; i < n_resamples; ++i) {
    //   sim = make_SpotSim(gtc.first, STdata.cb, fr, STdata.correction_table, 999 + i, n_forks);
    //   corrected_counts_sim.row(i) = NumericVector(sim.corrected_counts.begin(), sim.corrected_counts.end());
    //   CRPPV = compute_CRPPV(sim);
    //   CR_sim.row(i) = NumericVector(CRPPV.first.begin(), CRPPV.first.end());
    //   PPV_sim.row(i) = NumericVector(CRPPV.second.begin(), CRPPV.second.end());
    //   nll_sim[i] = compute_sim_nll(N_cells, STdata.bc_rates, gtc.second, sim.corrected_counts);
    //   if ((i + 1) % (n_resamples/10) == 0) {
    //     Rcpp::Rcout << "Simulation " << i + 1 << "/" << n_resamples << " complete." << std::endl;
    //   }
    // }
    
    return List::create(
      _["counts_provided"] = bc_counts,
      _["expected_counts"] = ecc,
      //_["counts_simulated"] = counts_simulated,
      //_["corrected_counts_sim"] = corrected_counts_sim,
      //_["CR_sim"] = CR_sim,
      //_["PPV_sim"] = PPV_sim,
      //_["nll_sim"] = nll_sim,
      _["rate10"] = fr.rate10,
      _["rate01"] = fr.rate01,
      _["corr"] = fr.corr//,
      //_["eval_history"] = STdata.eval_history
    );
  }
