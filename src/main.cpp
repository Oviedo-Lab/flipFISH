
// Rcpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <unordered_set>
#include <map>          // std::map, used by draw_structured_corr() to group bits by round/channel
#include <numeric>      // std::accumulate
#include "fork_utils.h" // run_forked_batches() -- shared fork/pipe map-reduce, Windows-guarded
#include <nlopt.hpp>    // L-BFGS and other gradient-based optimizers
#include <boost/math/distributions/normal.hpp>     // norm_cdf() -- native (no R callback) normal CDF/quantile
#include <boost/math/special_functions/owens_t.hpp> // bvn_cdf_uppertail() -- native bivariate normal CDF
#include <boost/math/distributions/poisson.hpp>     // dg_ppois_lower() -- native (no R callback) Poisson CDF
#include <boost/random/mersenne_twister.hpp>        // boost::random::mt19937 -- portable RNG engine
#include <boost/random/normal_distribution.hpp>     // boost::random::normal_distribution
#include <boost/random/bernoulli_distribution.hpp>  // boost::random::bernoulli_distribution
#include <boost/random/gamma_distribution.hpp>       // boost::random::gamma_distribution
#include <boost/random/poisson_distribution.hpp>     // boost::random::poisson_distribution
#include <boost/random/uniform_01.hpp>               // boost::random::uniform_01
#include <boost/random/uniform_int_distribution.hpp> // boost::random::uniform_int_distribution
using namespace Rcpp;
using namespace Eigen;

// Hamming distance
// int dist = __builtin_popcountll(a ^ b);

// Define aliases for compact code
using V64  = std::vector<uint64_t>;
using Vint = std::vector<int>;
using Vdbl = std::vector<double>;
using Vboo = std::vector<bool>;
using Vstr = std::vector<std::string>;
using LUT  = std::unordered_map<uint64_t, int>;
using LUTr = std::unordered_map<int, V64>;

/*
 * *********************************************************************************************************************
 * Basic data structures
 */

// Barcode codebook
struct Codebook {
    // ... Note: Experimental structure (rounds and channels) captured in FlipRatePriors
    int  N_bits;
    V64  barcodes; // list of barcodes
    Vstr species;  // list of species (gene or blank) names, order matches barcodes
    Vint blanks;   // list of indexes of barcodes and species, giving which are blanks
    Vint genes;    // same, but giving indexes of genes
  };

// Optimizer results tracking
struct EvalResults {
    Vdbl msle_hist;   // covering each optimizer step
    Vdbl ehc;         // expected hit count, for each barcode (best optimizer step)
    Vdbl ecc;         // expected corrected count, ... 
    Vdbl erc;         // expected read count, ... 
    int  n_evals = 0; // length of msle_hist
  };

// Flip rates
struct FlipRates {
    Vdbl rate10;      // P(1 -> 0) for each bit, so length = N_bits
    Vdbl rate01;      // P(0 -> 1) for each bit, so length = N_bits
    Vdbl log_rate10;  // Pre-computed log(rate10): log P(1->0) for each bit, used in TR
    Vdbl log_rate01;  // Pre-computed log(rate01): log P(0->1) for each bit, used in TR
    Vdbl corr1;       // Strict lower triangle of luminance noise correlations when bit i (i.e., column bit) is 1 
    Vdbl corr0;       // ... when bit i is 0
  }; 

// Flip rate priors
struct FlipRatePriors {
    double expected10           = 0.01; // expected rate at which 1 bits flip to 0s
    double expected01           = 0.05; // expected rate at which 0 bits flip to 1s
    double sd10                 = 0.1;  // standard deviation of the 1-to-0 flip rate
    double sd01                 = 0.2;  // standard deviation of the 0-to-1 flip rate
    double expectedcorr         = 0.0;  // expected baseline flip-correlation between bits, before considering channels and rounds
    double sdcorr               = 0.2;  // standard deviation of flip-correlation between bits, before considering channels and rounds
    double round_corr_sd        = 0.2;  // standard deviation fed into each round's own local correlation draw
    double channel_rho          = 0.5;  // equicorrelation value shared by every pair within a channel, before w_channel scaling (must be >= -1/(m-1) for the smallest channel group of size m)
    double w_round              = 0.7;  // fraction of each bit's variance attributed to round membership (0 disables the round tier)
    double w_channel            = 0.1;  // fraction attributed to channel membership (0 disables the channel tier); w_round + w_channel must be <= 1
    bool   use_latent_structure = false;
    Vdbl   sdcorr_pairs;
    Vint   latent_bit_round;
    Vint   latent_bit_channel;
    /*
     * Structured priors: 
     *    sdcorr_pairs: Optional per-pair override of sdcorr (same strict-lower-triangle pair index as corr1 and corr0,
     *      i.e. pair (i,j), j<i, at lower_tri_index(i, j)). Empty (default) means "not set" and falls back to sdcorr.
     *      Lets a caller express a tiered correlation prior (e.g. a wide prior for bit-pairs sharing an imaging round,
     *      a narrower one for pairs only sharing a color channel, and a tight near-zero prior for everything else)
     *      via R's build_corr_priors(), instead of one global sd for all pairs.
     *      
     *    latent_bit_round and latent_bit_channel: Two-level latent-factor structured-truth generative mode, used only 
     *      by draw_random_fr() -- NOT used by the fitting penalty above (that's sdcorr or sdcorr_pairs). Active 
     *      when use_latent_structure is true. latent_bit_round and latent_bit_channel are both length N_bits each, 
     *      arbitrary integer group codes; see R's build_structured_corrs()). Lets sim.benchmark() simulate a "true" 
     *      correlation structure organized by, e.g., imaging round and color channel, instead of i.i.d. draws across 
     *      all pairs.
     */
  };

// Spatial transcriptomics data
struct ST_data {
    // ST data, vectors of length N_barcodes
    Vdbl           bc_rates;            // Expected count, per cell
    Vdbl           bc_variance;         // Expected variance in count, among cells
    Vint           bc_counts;           // Total counts for each barcode across all cells
    // Codebook and look-up tables
    Codebook       cb;
    int            maxHam;              // Max correctable Hamming distance
    LUT            correction_table;
    LUTr           correction_table_inverted;
    // Parameter estimation
    EvalResults    eval_results;
    Vint           bc_counts_true;      // Expected true counts for each barcode, based on estimated misread loss
    double         best_msle;
    int            n_forks;             // For computing expected counts in parallel
    int            max_flips;           // For approximating expected counts estimate by ignoring highly unlikely misreads
    int            report_freq;
    // Regularization
    Vdbl           msle_weights;
    FlipRatePriors fr_priors;
    double         prior_weight       = 0.0;
    double         dispersion_weight  = 0.0;
    // Exact-match read counts
    bool           has_erc_obs        = false;
    Vint           bc_counts_erc_obs;
    double         erc_weight         = 0.0;
    // Feasibility penalty (shared with mQC_DG's MCMCSA objective; see dg_feasibility_penalty())
    double         feasibility_weight = 0.0;
    int            max_fs_bc          = 40;
    /*
     * Regularization terms:
     *
     *   1. per-barcode weights for the corrected-count MSLE term (used to up-weight blank barcodes,
     *   whose true hit count is known to be zero and which therefore carry information
     *   about misassignment/cross-talk that would otherwise be swamped by well-expressed genes).
     *
     *   2. a Gaussian prior pulling flip rates/correlations toward plausible values (mainly to regularize
     *   the otherwise unidentifiable pairwise bit-flip correlation terms).
     *
     *   3. dispersion_weight penalizes the empirical variance of corr1 (and, separately, corr0) across
     *   bit-pairs for deviating from sdcorr^2. Without this, a degenerate solution where all corr1[i] (or
     *   corr0[i]) collapse to the same value -- including all zero -- incurs no extra cost beyond the
     *   pointwise term, even though a nonzero sdcorr says we expect heterogeneity across bit-pairs.
     *
     *   4. feasibility_weight penalizes corr1/corr0 for requiring a large correction before they could
     *   be realized as an actual (PSD) per-barcode Gaussian correlation structure -- see
     *   dg_feasibility_penalty() -- so this fit stays comparable to what the DG method (mQC_DG) could
     *   actually simulate from the same corr1/corr0.
     *
     * Exact-match read counts:
     *
     *   Optional joint fit on exact-match (uncorrected) read counts, in addition to the corrected-count fit.
     *   This gives the optimizer a second, independent observable that constrains how much of each barcode's
     *   corrected count came from exact matches vs. correction, which is otherwise invisible to a
     *   corrected-count-only objective.
     */
  };

// Weights controlling the fit objective's up-weighting and regularization terms
struct Weights {
    double blank       = 1.0; // upweight blank barcodes' corrected/read-count MSLE contribution
    double prior       = 0.0; // Gaussian prior-penalty strength on flip rates/correlations (0 disables)
    double dispersion  = 0.0; // dispersion-penalty strength matching empirical variance of corr1/corr0 to sdcorr^2 (0 disables)
    double erc         = 0.0; // joint exact-match (uncorrected) read-count MSLE term strength (0 disables; ignored if bc_counts has no "counts_uncorrected" column)
    double feasibility = 0.0; // feasibility-penalty strength for corr1/corr0 needing a large correction to become realizable (0 disables)
    double max_blank   = 1.0; // mQC_DG() only: blank weight at the start of the FR_search/burn-in schedule
    double min_blank   = 1.0; // mQC_DG() only: blank weight by the end of burn-in (the schedule's held value)
    /*
     * Note: max_blank and min_blank are mQC_DG()-only (the blank-weight schedule's start and end); R always resolves 
     * these from the data before calling, since a sensible value depends on the gene/blank counts, but they still 
     * need some literal fallback here for a direct C++-side call that omits them.
     */
  };

// Spot simulation
struct SpotSim {
    Vint read_counts;      // For each barcode, the total number of spots read as that barcode,
    Vint corrected_counts; // ... corrected to that barcode
    Vint hit_counts;       // ... and correctly read as that barcode
  };

// One MCMCSA step's simulation-fit result
struct DGStepResult {
    double sim_mse;
    Vdbl   PPV_genes;        // length N_genes, may contain NaN
    Vdbl   read_counts;      // length N_barcodes, natural order (genes then blanks)
    Vdbl   corrected_counts; // length N_barcodes, natural order
    Vdbl   hit_counts;       // length N_barcodes, natural order
    /*
     * This is the scalar sim-mse objective that drives MCMCSA acceptance (and is reported out as misread.qc.DG()'s 
     * "msle"), this step's per-gene PPV (reported out as "PPV_by_step"), and this step's per-barcode read, corrected, 
     * and hit counts in natural STdata row order (genes and blanks, consumed by mQC_DG() to build "erctc_plus" -- never 
     * themselves returned to R, so the only things crossing the R/C++ boundary per step stay small). PPV_genes is NaN 
     * for any gene with zero corrected spots this step, which is expected and legitimate at small N_cells.
     */
  };

/*
 * *********************************************************************************************************************
 * Helper functions
 */

// Objective: mean squared log error
double compute_msle(
    const Vint& obs_counts,   // observed corrected counts
    const Vdbl& pred_counts,  // predicted corrected counts
    const Vdbl& weights       // per-barcode weight (e.g., to up-weight blanks)
  ) {
    int N_barcodes = obs_counts.size();
    if (pred_counts.size() != N_barcodes) {Rcpp::stop("obs_counts and pred_counts must be the same length.");}
    if (weights.size()     != N_barcodes) {Rcpp::stop("obs_counts and weights must be the same length.");}
    double msle       = 0.0;
    double weight_sum = 0.0;
    for (int b = 0; b < N_barcodes; ++b) {
      double le = std::log((pred_counts[b] + 1.0) / (static_cast<double>(obs_counts[b]) + 1.0)); // log(A) - log(B) = log(A/B)
      msle       += weights[b] * le * le;
      weight_sum += weights[b];
    }
    msle /= weight_sum;
    return msle;
  }

// String search returning indexes
Vint grep_idx(
    CharacterVector x,
    std::string     s,
    bool            neg = false
  ) {
    Vint idx;
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
    return idx;
  }

// Index of pair (i, j), j < i, within the strict lower triangle's canonical linear order.
inline int lower_tri_index(
    int i,
    int j
  ) {
    return i * (i - 1) / 2 + j;
   
    /*
     * This is the canonical linear order corr1/corr0/sdcorr_pairs all use throughout this file
     * (and that lower_tri_pairs_raw()/its R-facing export below enumerate in that same order).
     * Single source of truth for this formula: every other site that needs it (TR(),
     * dg_raw_gaussian_corr(), dg_feasibility_penalty(), draw_structured_corr()) calls this instead
     * of writing i*(i-1)/2+j itself, so the formula can't drift between call sites.
     */
  }

// Enumerate the strict lower triangle's canonical pair order for N_bits bits.
std::pair<Vint, Vint> lower_tri_pairs_raw(
    int N_bits
  ) {
    int  corr_free = N_bits * (N_bits - 1) / 2;
    Vint is(corr_free), js(corr_free);
    int  k = 0;
    for (int i = 1; i < N_bits; ++i) {
      for (int j = 0; j < i; ++j) {
        is[k] = i;
        js[k] = j;
        ++k;
      }
    }
    return {is, js};
   
    /*
     * Pair (i, j) sits at linear index lower_tri_index(i, j), and walking i = 1..N_bits-1,
     * j = 0..i-1 visits those indices in increasing order, so a plain running counter matches
     * lower_tri_index() without computing it. Returns that enumeration as two parallel 0-based
     * index vectors (is[k], js[k] give the pair at linear index k) of length N_bits*(N_bits-1)/2.
     * This is the one place that order is built; draw_random_fr() and correlation_from_params()
     * (both below) call it instead of re-deriving the same nested loop, as does lower_tri_pairs()
     * (the R-facing export just below this), so none of them can drift out of sync with each other
     * or with lower_tri_index().
     */
  }

// R-facing export of lower_tri_pairs_raw(), as a 1-based i/j index matrix.
// [[Rcpp::export]]
IntegerMatrix lower_tri_pairs(
    int N_bits
  ) {
    auto pairs     = lower_tri_pairs_raw(N_bits);
    int  corr_free = pairs.first.size();
    IntegerMatrix out(corr_free, 2);
    colnames(out) = CharacterVector::create("i", "j");
    for (int k = 0; k < corr_free; ++k) {
      out(k, 0) = pairs.first[k]  + 1;
      out(k, 1) = pairs.second[k] + 1;
    }
    return out;
   
    /*
     * The strict lower triangle's canonical pair order, as a corr_free x 2 integer matrix with
     * 1-based (R-style) "i"/"j" columns -- row k (1-based) gives the pair at linear index k,
     * matching corr1/corr0's/sdcorr_pairs' own layout. Lets build_corr_priors() (R/main.R)
     * enumerate these pairs directly instead of hand-rolling the matching nested loop.
     */
  }

// Convert flip rate vector into structured FR object
FlipRates pack_fr(
    const Vdbl& rates_plus_corr,
    int         N_bits
  ) {
    // Initialize flip-rate structure
    FlipRates fr;
    // Assign flip rates
    fr.rate10.assign(rates_plus_corr.begin(), rates_plus_corr.begin() + N_bits);
    fr.rate01.assign(rates_plus_corr.begin() + N_bits, rates_plus_corr.begin() + 2*N_bits);
    // Find number of correlations (strict lower triangle) ... corr(i,j) = corrXXXX(lower_tri_index(i, j))
    int N_corrs = N_bits * (N_bits - 1) / 2;
    // Assign correlations
    fr.corr1.assign(rates_plus_corr.begin() + 2*N_bits, rates_plus_corr.begin() + 2*N_bits + N_corrs);
    fr.corr0.assign(rates_plus_corr.begin() + 2*N_bits + N_corrs, rates_plus_corr.end());
    // Precompute log flip rates used in TR
    fr.log_rate10 = Vdbl(N_bits, 0.0);
    fr.log_rate01 = Vdbl(N_bits, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      fr.log_rate10[i] = std::log(fr.rate10[i]);
      fr.log_rate01[i] = std::log(fr.rate01[i]);
    }
    return fr;
  }

// Extract flip rate priors from input
FlipRatePriors pack_fr_priors(
    const List& fliprate_priors
  ) {
    FlipRatePriors fr_priors; // Initializes with all defaults
    if (fliprate_priors.containsElementNamed("expected10"))
      fr_priors.expected10           = as<double>(fliprate_priors["expected10"]);
    if (fliprate_priors.containsElementNamed("expected01"))
      fr_priors.expected01           = as<double>(fliprate_priors["expected01"]);
    if (fliprate_priors.containsElementNamed("sd10"))
      fr_priors.sd10                 = as<double>(fliprate_priors["sd10"]);
    if (fliprate_priors.containsElementNamed("sd01"))
      fr_priors.sd01                 = as<double>(fliprate_priors["sd01"]);
    if (fliprate_priors.containsElementNamed("expectedcorr"))
      fr_priors.expectedcorr         = as<double>(fliprate_priors["expectedcorr"]);
    if (fliprate_priors.containsElementNamed("sdcorr"))
      fr_priors.sdcorr               = as<double>(fliprate_priors["sdcorr"]);
    // Optional tiered prior: per-pair sdcorr override (see FlipRatePriors::sdcorr_pairs)
    if (fliprate_priors.containsElementNamed("sdcorr_pairs")) {
      NumericVector sp               = as<NumericVector>(fliprate_priors["sdcorr_pairs"]);
      fr_priors.sdcorr_pairs         = Vdbl(sp.begin(), sp.end());
    }
    // Optional structured-truth generative mode (see FlipRatePriors::latent_bit_round etc.)
    if (fliprate_priors.containsElementNamed("latent_bit_round")) {
      IntegerVector lr               = as<IntegerVector>(fliprate_priors["latent_bit_round"]);
      fr_priors.latent_bit_round     = Vint(lr.begin(), lr.end());
    }
    if (fliprate_priors.containsElementNamed("latent_bit_channel")) {
      IntegerVector lc               = as<IntegerVector>(fliprate_priors["latent_bit_channel"]);
      fr_priors.latent_bit_channel   = Vint(lc.begin(), lc.end());
    }
    if (fliprate_priors.containsElementNamed("round_corr_sd"))
      fr_priors.round_corr_sd        = as<double>(fliprate_priors["round_corr_sd"]);
    if (fliprate_priors.containsElementNamed("channel_rho"))
      fr_priors.channel_rho          = as<double>(fliprate_priors["channel_rho"]);
    if (fliprate_priors.containsElementNamed("w_round"))
      fr_priors.w_round              = as<double>(fliprate_priors["w_round"]);
    if (fliprate_priors.containsElementNamed("w_channel"))
      fr_priors.w_channel            = as<double>(fliprate_priors["w_channel"]);
    if (fliprate_priors.containsElementNamed("use_latent_structure"))
      fr_priors.use_latent_structure = as<bool>(fliprate_priors["use_latent_structure"]);
    return fr_priors;
  }

// Extract weights from input (same containsElementNamed pattern as pack_fr_priors())
Weights pack_weights(
    const List& weights
  ) {
    Weights w; // Initializes with all defaults
    if (weights.containsElementNamed("blank"))       w.blank       = as<double>(weights["blank"]);
    if (weights.containsElementNamed("prior"))       w.prior       = as<double>(weights["prior"]);
    if (weights.containsElementNamed("dispersion"))  w.dispersion  = as<double>(weights["dispersion"]);
    if (weights.containsElementNamed("erc"))         w.erc         = as<double>(weights["erc"]);
    if (weights.containsElementNamed("feasibility")) w.feasibility = as<double>(weights["feasibility"]);
    if (weights.containsElementNamed("max_blank"))   w.max_blank   = as<double>(weights["max_blank"]);
    if (weights.containsElementNamed("min_blank"))   w.min_blank   = as<double>(weights["min_blank"]);
    return w;
  }

// Invert lookup table
LUTr invert_lookup_table(
    const LUT& lut
  ) {
    LUTr inverted;
    for (const auto& kv : lut) {
      uint64_t barcode = kv.first;
      int      label   = kv.second;
      inverted[label].push_back(barcode);
    }
    return inverted;
  }

// Pack vector barcode into integer
uint64_t pack(
    Vint bits
  ) {
    uint64_t packed = 0;
    for (int i = 0; i < bits.size(); ++i) {packed |= (uint64_t(bits[i]) << i);} 
    return packed;
  }

// Helper for generating barcodes near x
void generate_neighbors(
    uint64_t x,
    int      n_bits,
    int      dist,
    int      start_bit,
    V64&     out
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

// Generate all barcodes near x
V64 neighbors(
    uint64_t x,
    int      n_bits,
    int      hamming_dist
  ) {
    V64 out;
    generate_neighbors(
      x,
      n_bits,
      hamming_dist,
      0, // start bit
      out
    );
    return out;
  }

// Create valid correlation matrix from random parameter draw
Eigen::MatrixXd correlation_from_params(
    const Vdbl& theta,
    int         n
  ) {
    Eigen::MatrixXd L = Eigen::MatrixXd::Identity(n, n);
    auto pairs = lower_tri_pairs_raw(n); // same canonical (i, j) order in which theta is laid out
    for (size_t k = 0; k < pairs.first.size(); ++k) {
      L(pairs.first[k], pairs.second[k]) = theta[k];
    }
    Eigen::MatrixXd S  = L * L.transpose();
    Eigen::VectorXd sd = S.diagonal().cwiseSqrt();
    Eigen::MatrixXd R  = S.array().colwise() / sd.array();
    R                  = R.array().rowwise() / sd.transpose().array();
    return R;
    // ... written by ChatGPT
  }

/*
 * *********************************************************************************************************************
 * Dichotomized-Gaussian (DG) correlation-matrix helpers
 *
 * corr1 and corr0 (FlipRates) are correlations of the *dichotomized* (post-threshold, 0/1) bit-read
 * outcomes -- the same quantity TR() already treats them as. Simulating from an actual dichotomized
 * Gaussian needs the underlying *Gaussian* correlation that reproduces a target dichotomized
 * correlation once thresholded; that inversion is what dg_find_sigma_RootBisection() does, via the
 * bivariate-normal-CDF identity Cov(1[X>h],1[Y>k]) = P(X>h,Y>k) - P(X>h)P(Y>k).
 *
 * bvn_cdf_uppertail() computes P(X>h, Y>k) natively (Owen's T function, boost::math::owens_t) rather
 * than calling back into R's mvtnorm::pmvnorm -- the bivariate case has a closed form, unlike the
 * general N-dimensional case mvtnorm handles, so this stays fully native and fast inside the hot
 * per-pair inversion loop. Verified against R's mvtnorm::pmvnorm (exact agreement, |diff| < 1e-10)
 * over a grid spanning both signs of h and k and rho in [-0.9, 0.9].
 */

static boost::math::normal_distribution<double> dg_std_normal(0.0, 1.0);

// Normal CDF, with inverse (native; no R callback)
double norm_cdf(
    double x,
    double mu,
    double sd,
    bool   inverse
  ) {
    if (inverse) {
      double xc  = x;
      double eps = std::numeric_limits<double>::epsilon();
      if (xc < eps)       {xc = eps;}
      if (xc > 1.0 - eps) {xc = 1.0 - eps;}
      boost::math::normal_distribution<double> dist(mu, sd);
      return boost::math::quantile(dist, xc);
    } else {
      boost::math::normal_distribution<double> dist(mu, sd);
      return boost::math::cdf(dist, x);
    }
  }

// Bivariate standard-normal upper-tail probability P(X>h, Y>k) with correlation rho.
double bvn_cdf_uppertail(
    double h,
    double k,
    double rho
  ) {
    double denom      = std::sqrt(1.0 - rho * rho);
    double a          = (k - rho * h) / (h * denom);
    double b          = (h - rho * k) / (k * denom);
    double delta      = h * k < 0.0 || (h * k == 0.0 && (h + k) < 0.0) ? 0.5 : 0.0;
    double Phi_h      = boost::math::cdf(dg_std_normal, h);
    double Phi_k      = boost::math::cdf(dg_std_normal, k);
    double phi2_lower = 0.5 * (Phi_h + Phi_k) - boost::math::owens_t(h, a) - boost::math::owens_t(k, b) - delta;
    return 1.0 - Phi_h - Phi_k + phi2_lower;
   
    /*
     * Via the Owen (1956) lower-tail identity Phi2(h,k;rho) = 0.5*(Phi(h)+Phi(k)) - T(h,a) - T(k,b)
     * - delta, then upper = 1 - Phi(h) - Phi(k) + Phi2(h,k;rho). Assumes h != 0 and k != 0 (always
     * true here: h/k are qnorm() of a flip rate, which is never exactly 0.5).
     */
  }

// Residual for dg_find_sigma_RootBisection(): desired_cov minus the achieved dichotomized covariance.
double dg_sigma_formula_scalar(
    double threshold_i,
    double threshold_j,
    double desired_cov,
    double rho
  ) {
    double Phi2_upper = bvn_cdf_uppertail(threshold_i, threshold_j, rho);
    double p_i        = 1.0 - boost::math::cdf(dg_std_normal, threshold_i); // P(read_i = 1)
    double p_j        = 1.0 - boost::math::cdf(dg_std_normal, threshold_j); // P(read_j = 1)
    return desired_cov - Phi2_upper + p_i * p_j;
  }

// Find the Gaussian correlation rho whose dichotomization reproduces a desired covariance.
double dg_find_sigma_RootBisection(
    double threshold_i,
    double threshold_j,
    double desired_cov
  ) {
    const int    max_iter = 60;
    const double tol      = 1e-10;
    double       lo       = -0.999;
    double       hi       = 0.999;
    double       f_lo     = dg_sigma_formula_scalar(threshold_i, threshold_j, desired_cov, lo);
    double       f_hi     = dg_sigma_formula_scalar(threshold_i, threshold_j, desired_cov, hi);
    if (std::fabs(f_lo) < tol) {return lo;}
    if (std::fabs(f_hi) < tol) {return hi;}
    if (f_lo * f_hi > 0.0)     {return 0.0;} // target covariance infeasible for these marginals
    double mid   = 0.0;
    double f_mid = std::numeric_limits<double>::infinity();
    int    iter  = 0;
    while (std::fabs(f_mid) > tol && iter < max_iter) {
      mid   = 0.5 * (lo + hi);
      f_mid = dg_sigma_formula_scalar(threshold_i, threshold_j, desired_cov, mid);
      if ((f_mid > 0.0) == (f_lo > 0.0)) {lo = mid; f_lo = f_mid;} else {hi = mid;}
      iter++;
    }
    return mid;
   
    /*
     * desired_cov is the raw covariance of the two dichotomized (0/1) outcomes at thresholds
     * (threshold_i, threshold_j), i.e. desired_corr * sqrt(p_i*(1-p_i)) * sqrt(p_j*(1-p_j)). Falls
     * back to rho = 0 when desired_cov lies outside what's achievable for these marginals at all
     * (the classic Frechet-Hoeffding bound on correlation given fixed marginals) -- same convention
     * as the historical implementation this is restored from. Verified end-to-end against the exact
     * achieved-vs-target correlation (not just the residual formula) for asymmetric marginals
     * (|diff| < 1e-10 whenever the target was feasible).
     */
  }

// Per-bit threshold tau such that "read = 1" iff Z > tau (Z ~ N(0,1)), given this bit's true value.
double dg_bit_threshold(
    int    bit_value,
    double rate10_i,
    double rate01_i
  ) {
    return bit_value == 1 ? norm_cdf(rate10_i, 0.0, 1.0, true) : norm_cdf(1.0 - rate01_i, 0.0, 1.0, true);
    // Matches TR()'s own flip-rate semantics: P(read=0|true=1) = rate10 and P(read=1|true=0) = rate01.
  }

// Build the raw (not yet guaranteed valid/PSD) Gaussian correlation matrix for one barcode's bits.
Eigen::MatrixXd dg_raw_gaussian_corr(
    const Vint&      barcode_bits, // length N_bits, values 0/1
    const FlipRates& fr,
    const Vdbl&      thresholds    // length N_bits, this barcode's own dg_bit_threshold() per bit
  ) {
    int             N_bits = barcode_bits.size();
    Eigen::MatrixXd Sigma  = Eigen::MatrixXd::Identity(N_bits, N_bits);
    for (int i = 1; i < N_bits; ++i) {
      for (int j = 0; j < i; ++j) {
        int    k            = lower_tri_index(i, j);
        double desired_corr = barcode_bits[i] == 1 ? fr.corr1[k] : fr.corr0[k];
        double p_i          = 1.0 - boost::math::cdf(dg_std_normal, thresholds[i]);
        double p_j          = 1.0 - boost::math::cdf(dg_std_normal, thresholds[j]);
        double desired_cov  = desired_corr * std::sqrt(p_i * (1.0 - p_i)) * std::sqrt(p_j * (1.0 - p_j));
        double rho          = dg_find_sigma_RootBisection(thresholds[i], thresholds[j], desired_cov);
        Sigma(i, j) = rho;
        Sigma(j, i) = rho;
      }
    }
    return Sigma;
   
    /*
     * Selects corr1[ij]/corr0[ij] per pair using the same convention TR() already uses (keyed on
     * the larger-indexed bit's value in this barcode -- see TR(), "fr.corr1[ij] : fr.corr0[ij]").
     */
  }

// Native port of R's Matrix::nearPD() -- Higham (2002)'s alternating-projections algorithm.
Eigen::MatrixXd dg_nearPD_native(
    const Eigen::MatrixXd& x,
    bool                   corr      = true,
    bool                   keepDiag  = true,
    bool                   do2eigen  = true,
    bool                   doDykstra = true,
    int                    maxit     = 100,
    double                 eig_tol   = 1e-6,
    double                 conv_tol  = 1e-7,
    double                 posd_tol  = 1e-8
  ) {
    int  n        = x.rows();
    auto inf_norm = [](const Eigen::MatrixXd& M) {
      return M.cwiseAbs().rowwise().sum().maxCoeff(); // R's norm(., "I"): max absolute row sum
    };
   
    Eigen::VectorXd diagX0    = x.diagonal();         // only read back when keepDiag && !corr
    Eigen::MatrixXd D_S       = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd X         = x;
    int             iter      = 0;
    bool            converged = false;
    while (iter < maxit && !converged) {
      Eigen::MatrixXd Y         = X;
      Eigen::MatrixXd Rk        = doDykstra ? (Y - D_S) : Y;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Rk);
      Eigen::VectorXd d         = es.eigenvalues();   // ascending
      Eigen::MatrixXd Q         = es.eigenvectors();
      double          threshold = eig_tol * d(n - 1); // d(n-1) here == R's d[1] (largest)
      Eigen::MatrixXd Xnew      = Eigen::MatrixXd::Zero(n, n);
      bool            any_kept  = false;
      for (int k = 0; k < n; ++k) {
        if (d(k) > threshold) {
          any_kept = true;
          Xnew    += d(k) * (Q.col(k) * Q.col(k).transpose());
        }
      }
      if (!any_kept) {Rcpp::stop("dg_nearPD_native: matrix seems negative semi-definite.");}
      X = Xnew;
      if (doDykstra)     {D_S = X - Rk;}
      if (corr)          {X.diagonal().setOnes();}
      else if (keepDiag) {X.diagonal() = diagX0;}
      double conv = inf_norm(Y - X) / inf_norm(Y);
      ++iter;
      converged = (conv <= conv_tol);
    }
    if (do2eigen) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(X);
      Eigen::VectorXd d = es2.eigenvalues(); // ascending
      Eigen::MatrixXd Q = es2.eigenvectors();
      double Eps = posd_tol * std::fabs(d(n - 1)); // d(n-1) here == R's d[1] (largest)
      if (d(0) < Eps) { // d(0) here == R's d[n] (smallest)
        Eigen::VectorXd dclipped = d;
        for (int k = 0; k < n; ++k) {if (dclipped(k) < Eps) {dclipped(k) = Eps;}}
        Eigen::VectorXd o_diag = X.diagonal();
        X = Q * dclipped.asDiagonal() * Q.transpose();
        Eigen::VectorXd Dscale(n);
        for (int i = 0; i < n; ++i) {Dscale(i) = std::sqrt(std::max(Eps, o_diag(i)) / X(i, i));}
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {X(i, j) *= Dscale(i) * Dscale(j);}
        }
      }
      if (corr)          {X.diagonal().setOnes();}
      else if (keepDiag) {X.diagonal() = diagX0;}
    }
    return X;
   
    /*
     * Matches Matrix::nearPD()'s own defaults/logic exactly (eigen() there returns eigenvalues
     * descending; Eigen::SelfAdjointEigenSolver returns them ascending, so every reference to R's
     * d[1]/d[n] above is reindexed to this function's d(n-1)/d(0) accordingly).
     *
     * Replaces an R callback that used to sit in make_SpotSim_DG()'s hot loop -- called once per
     * distinct active barcode on every single MCMCSA proposal, accepted or not, so each
     * R-interpreter round trip was paid many times per step once the annealing's acceptance rate
     * collapsed (see draw_random_fr() et al. for the RNG side of that same loop). Being native also
     * removes the fork-unsafety the R callback used to impose on dg_precompute_barcode_noise() (see
     * make_SpotSim_DG()'s own comment on why that precompute pass used to have to run serially, in
     * the parent process, before forking).
     *
     * Verified against R's Matrix::nearPD(corr=TRUE, keepDiag=TRUE, do2eigen=TRUE) (its call
     * signature here) over random indefinite/invalid correlation-like matrices, including
     * dg_raw_gaussian_corr()'s actual output shape (up to 28x28): agrees to the same magnitude as
     * float round-off, with the same number of iterations and the same achieved correlations.
     */
  }

// Project a (possibly invalid) correlation matrix to the nearest valid one (see dg_nearPD_native()).
// ... corr/keepDiag/do2eigen = TRUE, as this package has always used
Eigen::MatrixXd dg_nearest_corr_matrix(const Eigen::MatrixXd& raw) {return dg_nearPD_native(raw, true, true, true);}

// Dichotomized correlation a single pair would achieve (ignoring nearPD's cross-pair redistribution).
double dg_pairwise_achieved_corr(
    double threshold_i,
    double threshold_j,
    double desired_corr,
    double sd_i,
    double sd_j
  ) {
    double desired_cov  = desired_corr * sd_i * sd_j;
    double rho          = dg_find_sigma_RootBisection(threshold_i, threshold_j, desired_cov);
    double p_i          = 1.0 - boost::math::cdf(dg_std_normal, threshold_i);
    double p_j          = 1.0 - boost::math::cdf(dg_std_normal, threshold_j);
    double achieved_cov = bvn_cdf_uppertail(threshold_i, threshold_j, rho) - p_i * p_j;
    return achieved_cov / (sd_i * sd_j);
   
    /*
     * Returns the raw pairwise map alone, for a given desired correlation. Used both to build each
     * barcode's raw matrix and, cheaply, to estimate dg_feasibility_penalty()'s gradient without
     * ever re-calling nearPD.
     */
  }

// How far corr1/corr0 would have to be distorted to become a valid (PSD) per-barcode Gaussian correlation structure
double dg_feasibility_penalty(
    const FlipRates& fr,
    const Codebook&  cb,
    int              max_barcodes = 40,
    bool             do_grad      = false,
    Vdbl*            grad         = nullptr
  ) {
    int N_bits     = fr.rate10.size();
    int corr_free  = N_bits * (N_bits - 1) / 2;
    if (corr_free == 0) {return 0.0;}
   
    int N_barcodes = cb.barcodes.size();
    int step       = std::max(1, N_barcodes / max_barcodes);
    
    Vdbl thresholds_1(N_bits), thresholds_0(N_bits); // this bit's threshold when it reads 1 vs 0
    for (int i = 0; i < N_bits; ++i) {
      thresholds_1[i] = dg_bit_threshold(1, fr.rate10[i], fr.rate01[i]);
      thresholds_0[i] = dg_bit_threshold(0, fr.rate10[i], fr.rate01[i]);
    }
   
    const double eps      = 1e-3;
    double       total_sq = 0.0;
    int          n_terms  = 0;
    for (int b = 0; b < N_barcodes; b += step) {
      Vint barcode_bits(N_bits);
      Vdbl thresholds(N_bits);
      for (int i = 0; i < N_bits; ++i) {
        barcode_bits[i] = (cb.barcodes[b] >> i) & 1ULL;
        thresholds[i]   = barcode_bits[i] == 1 ? thresholds_1[i] : thresholds_0[i];
      }
      Eigen::MatrixXd raw       = dg_raw_gaussian_corr(barcode_bits, fr, thresholds);
      Eigen::MatrixXd projected = dg_nearest_corr_matrix(raw);
      for (int i = 1; i < N_bits; ++i) {
        for (int j = 0; j < i; ++j) {
          int    k             = lower_tri_index(i, j);
          bool   bit_i_is_1    = barcode_bits[i] == 1;
          double desired_corr  = bit_i_is_1 ? fr.corr1[k] : fr.corr0[k];
          double p_i           = 1.0 - boost::math::cdf(dg_std_normal, thresholds[i]);
          double p_j           = 1.0 - boost::math::cdf(dg_std_normal, thresholds[j]);
          double sd_i          = std::sqrt(p_i * (1.0 - p_i));
          double sd_j          = std::sqrt(p_j * (1.0 - p_j));
          double achieved_cov  = bvn_cdf_uppertail(thresholds[i], thresholds[j], projected(i, j)) - p_i * p_j;
          double achieved_corr = achieved_cov / (sd_i * sd_j);
          double diff          = achieved_corr - desired_corr;
          total_sq            += diff * diff;
          n_terms++;
         
          if (do_grad) {
            double ach_plus     = dg_pairwise_achieved_corr(thresholds[i], thresholds[j], desired_corr + eps, sd_i, sd_j);
            double ach_minus    = dg_pairwise_achieved_corr(thresholds[i], thresholds[j], desired_corr - eps, sd_i, sd_j);
            double sensitivity  = (ach_plus - ach_minus) / (2.0 * eps); // d(achieved_unprojected)/d(desired)
            double dpenalty     = 2.0 * diff * (sensitivity - 1.0);     // d(diff^2)/d(desired)
            int    param_idx    = bit_i_is_1 ? (2*N_bits + k) : (2*N_bits + corr_free + k);
            (*grad)[param_idx] += dpenalty; // divided by n_terms once, after the loop
          }
        }
      }
    }
    if (n_terms == 0) {return 0.0;}
    if (do_grad) {
      for (int k = 0; k < 2*corr_free; ++k) {(*grad)[2*N_bits + k] /= static_cast<double>(n_terms);}
    }
    return total_sq / static_cast<double>(n_terms);
   
    /*
     * Measured in dichotomized-correlation units so it's directly comparable to corr1/corr0
     * themselves: for a deterministic subsample of up to max_barcodes distinct barcodes, build the
     * raw Gaussian correlation matrix, project it to the nearest valid one, convert each
     * off-diagonal entry's *achieved* dichotomized correlation back out (same forward formula used
     * to build it), and accumulate the squared difference from what corr1/corr0 actually asked
     * for. Shared by mQC_msle() (penalizes the analytic fit for drifting somewhere DG couldn't
     * realize) and mQC_DG()'s MCMCSA objective (added directly as part of the per-step
     * simulation's own matrix assembly).
     *
     * Gradient (when do_grad): exact wrt the nearPD projection would require differentiating
     * through an iterative alternating-projections algorithm, which has no simple closed form.
     * Instead, each pair's contribution to d(achieved)/d(desired) is estimated locally via
     * dg_pairwise_achieved_corr() alone -- i.e. ignoring nearPD's cross-pair redistribution --
     * which costs one extra cheap (no-nearPD, no-R-callback) root-bisection per pair instead of
     * re-running the whole expensive projection. This still pushes optimization away from grossly
     * infeasible per-pair targets (where the local map is flat/clamped and the approximation is
     * exact); it under-credits distortion that arises *only* from joint infeasibility across
     * pairs, which the nearPD-based value term still captures and reports.
     */
  }

/*
 * *********************************************************************************************************************
 * Functions for spot decoding
 */

// Extract unique Hamming distances between barcodes
Vint unique_Hamming(
    const V64& barcodes
  ) {
    int N_barcodes = barcodes.size();
    Vint hamming_distances;
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
    bool verbose = true
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
      Vint bits(N_bits);
      for (int j = 0; j < N_bits; ++j) {bits[j] = codebook(i, j);}
      cb.barcodes.push_back(pack(bits));
      cb.species.push_back(as<std::string>(species[i]));
    }
    // Find all unique Hamming distances between barcodes 
    Vint hamming_distances = unique_Hamming(cb.barcodes);
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
    Codebook cb                = pack_codebook(codebook, false);
    Vint     hamming_distances = unique_Hamming(cb.barcodes);
    return wrap(hamming_distances);
  }

// Build table used to correct read spot barcodes
LUT build_correction_table(
    const Codebook& cb,
    int maxHam
  ) {
    LUT correction_table;
    std::unordered_set<uint64_t> ambiguous;
    for (size_t i = 0; i < cb.barcodes.size(); ++i) {
      correction_table[cb.barcodes[i]] = i; // Exact match
      // Generate all barcodes within maxHam
      for (int d = 1; d <= maxHam; ++d) {
        V64 n = neighbors(cb.barcodes[i], cb.N_bits, d);
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
 * Functions to analytically (conditional probability) compute expected counts from flip rates and true counts
 */

// Draw one set of bit-flip correlations
static Vdbl draw_structured_corr(
    int                     N_bits,
    const Vint&             bit_round,
    const Vint&             bit_channel,
    double                  round_corr_sd,
    double                  channel_rho,
    double                  w_round,
    double                  w_channel,
    boost::random::mt19937& rng
  ) {
    int  corr_free = N_bits * (N_bits - 1) / 2;
    Vdbl corr(corr_free, 0.0);
   
    // Group bit indices (0-based) by round and by channel
    std::map<int, Vint> rounds, channels;
    for (int i = 0; i < N_bits; ++i) {
      rounds[bit_round[i]].push_back(i);
      channels[bit_channel[i]].push_back(i);
    }
   
    // Round tier
    boost::random::normal_distribution<double> theta_draw(0.0, round_corr_sd);
    for (auto& kv : rounds) {
      const Vint& members = kv.second;
      int         m       = (int)members.size();
      if (m < 2 || w_round <= 0.0) continue;
      Vdbl theta(m * (m - 1) / 2);
      for (double& t : theta) t = theta_draw(rng);
      auto local_mat = correlation_from_params(theta, m);
      for (int a = 1; a < m; ++a) {
        for (int b = 0; b < a; ++b) {
          int hi  = std::max(members[a], members[b]);
          int lo  = std::min(members[a], members[b]);
          int idx = lower_tri_index(hi, lo);
          corr[idx] = w_round * local_mat(a, b);
        }
      }
    }
   
    // Channel tier: one shared, homogeneous value per channel, skipping any pair that already got
    // a round-tier value above (for this package's intended use that never happens -- a round's
    // bits always differ in channel -- but the check keeps this safe to reuse on a codebook where
    // it might).
    for (auto& kv : channels) {
      const Vint& members = kv.second;
      int         m       = (int)members.size();
      if (m < 2 || w_channel <= 0.0) continue;
      for (int a = 1; a < m; ++a) {
        for (int b = 0; b < a; ++b) {
          if (bit_round[members[a]] == bit_round[members[b]]) continue;
          int hi  = std::max(members[a], members[b]);
          int lo  = std::min(members[a], members[b]);
          int idx = lower_tri_index(hi, lo);
          corr[idx] = w_channel * channel_rho;
        }
      }
    }
   
    return corr;
    
    /*
     * Draw one set of bit-flip correlations (corr1 or corr0, depending on caller) via the two-level 
     * latent-factor model on FlipRatePriors::latent_bit_round/latent_bit_channel. This is provably a 
     * valid correlation structure -- not approximately, not pending a PSD check -- because it is
     * literally the covariance of
     *    X_i = sqrt(w_round)*(round-local term) + sqrt(w_channel)*(channel-shared term)
     *        + sqrt(1 - w_round - w_channel)*(idiosyncratic term),
     * three independent, zero-mean, unit-variance pieces, so Var(X_i) = 1 always, and for i != j:
     *    - sharing a round:  Cov(X_i,X_j) = w_round * (that round's own local correlation for the pair)
     *    - sharing a channel (not a round): Cov(X_i,X_j) = w_channel * channel_rho
     *    - neither: Cov(X_i,X_j) = 0 exactly (no cross-tier leakage at all)
     * The round-local term reuses correlation_from_params() (the package's existing Cholesky trick)
     * but restricted to just that round's own bits, so it's fully flexible (can be positive or
     * negative, heterogeneous per pair) without touching anything outside that round. The channel term
     * is a single shared (homogeneous) equicorrelation value per channel, which is why channel_rho
     * must satisfy the standard equicorrelation PSD bound, channel_rho >= -1/(m-1), for the smallest
     * channel group of size m (not checked here; pick a non-negative value and this never binds).
     * 
     * Initial draft of this function written by Claude Sonnet 5. 
     */
  }

// Select random flip rates and bit-flip correlations
Vdbl draw_random_fr(
    const ST_data&        STdata,
    const FlipRatePriors& fr_prior,
    unsigned int          seed = 12345 // Fixed default preserves prior deterministic behavior at existing call sites
  ) {
    // Initialize normal distributions
    boost::random::mt19937                     rng(seed);
    boost::random::normal_distribution<double> norm10(std::sqrt(fr_prior.expected10), fr_prior.sd10);
    boost::random::normal_distribution<double> norm01(std::sqrt(fr_prior.expected01), fr_prior.sd01);
    boost::random::normal_distribution<double> normCorr(fr_prior.expectedcorr, fr_prior.sdcorr);
    
    // Randomly select flip rates
    int    N_bits     = STdata.cb.N_bits;
    int    corr_free  = N_bits * (N_bits - 1) / 2;
    size_t n          = 2*N_bits + 2*corr_free;
    double near_one   = 1.0 - std::numeric_limits<double>::epsilon();
    double FRsqrt;
    Vdbl FR(n, 0.0);
    for (int i = 0; i < N_bits; ++i) {
      FRsqrt         = norm10(rng);
      FR[i]          = std::min(FRsqrt*FRsqrt, near_one);
      FRsqrt         = norm01(rng);
      FR[N_bits + i] = std::min(FRsqrt*FRsqrt, near_one);
    }
    // ... and random bit-flip correlations, by one of two routes:
    Vdbl corr1_final(corr_free), corr0_final(corr_free);
    if (fr_prior.use_latent_structure) {
      // Structured-truth generative mode (see FlipRatePriors::latent_bit_round and
      // draw_structured_corr() above). Draws corr1/corr0 directly in final, already-valid form --
      // no global Cholesky-normalize step, so no cross-tier leakage.
      if ((int)fr_prior.latent_bit_round.size()   != N_bits ||
          (int)fr_prior.latent_bit_channel.size() != N_bits) {
        Rcpp::stop("fr_prior.latent_bit_round/latent_bit_channel must each have length N_bits.");
      }
      if (fr_prior.w_round + fr_prior.w_channel > 1.0) {
        Rcpp::stop("fr_prior.w_round + fr_prior.w_channel must be <= 1.");
      }
      corr1_final = draw_structured_corr(
        N_bits, fr_prior.latent_bit_round, fr_prior.latent_bit_channel,
        fr_prior.round_corr_sd, fr_prior.channel_rho, fr_prior.w_round, fr_prior.w_channel, rng
      );
      corr0_final = draw_structured_corr(
        N_bits, fr_prior.latent_bit_round, fr_prior.latent_bit_channel,
        fr_prior.round_corr_sd, fr_prior.channel_rho, fr_prior.w_round, fr_prior.w_channel, rng
      );
    } else {
      // Original route: draw each pair (i.i.d., or per-pair sd when fr_prior.sdcorr_pairs is set --
      // see R's build_corr_priors()), then run the result through the package's global Cholesky
      // trick (correlation_from_params) to guarantee a valid correlation matrix. normCorr's
      // distribution object is reused via its param_type overload so per-pair sd costs nothing
      // extra when sdcorr_pairs is empty (the common case).
      bool has_pair_sd = !fr_prior.sdcorr_pairs.empty();
      if (has_pair_sd && (int)fr_prior.sdcorr_pairs.size() != corr_free) {
        Rcpp::stop("fr_prior.sdcorr_pairs must have length N_bits*(N_bits-1)/2 (one sd per bit pair).");
      }
      Vdbl corr1(corr_free);
      Vdbl corr0(corr_free);
      for (int i = 0; i < corr_free; ++i) {
        double draw1, draw0;
        if (has_pair_sd) {
          boost::random::normal_distribution<double>::param_type p(fr_prior.expectedcorr, fr_prior.sdcorr_pairs[i]);
          draw1 = normCorr(rng, p);
          draw0 = normCorr(rng, p);
        } else {
          draw1 = normCorr(rng);
          draw0 = normCorr(rng);
        }
        corr1[i] = std::max(std::min(draw1, near_one), -near_one);
        corr0[i] = std::max(std::min(draw0, near_one), -near_one);
      }
      auto corr1_mat = correlation_from_params(corr1, N_bits);
      auto corr0_mat = correlation_from_params(corr0, N_bits);
      auto pairs     = lower_tri_pairs_raw(N_bits);
      for (size_t m = 0; m < pairs.first.size(); ++m) {
        corr1_final[m] = corr1_mat(pairs.first[m], pairs.second[m]);
        corr0_final[m] = corr0_mat(pairs.first[m], pairs.second[m]);
      }
    }
    // ... fill into FR (both routes produce corr1_final/corr0_final already in lower_tri_index()'s
    // canonical order, matching the loop order below)
    int k = N_bits * 2;
    for (int idx = 0; idx < corr_free; ++idx) {FR[k++] = corr1_final[idx];}
    for (int idx = 0; idx < corr_free; ++idx) {FR[k++] = corr0_final[idx];}
    return FR;
  }

// Parse count data and codebook into ST_data structure
ST_data load_STdata(
    NumericMatrix bc_count_data,
    IntegerMatrix codebook,
    int           max_flips,
    int           maxHam,
    bool          verbose
  ) {
    
    // Load codebook as packed integers
    Codebook cb = pack_codebook(codebook, verbose);
    
    // Check max_flips value: < 1 (e.g. an explicit 0, or a direct Rcpp-level call bypassing the R
    // wrapper's own default) means "no limit", i.e. every misread the correction table admits gets
    // a full TR() evaluation. This is almost always prohibitively slow on a real codebook -- cost
    // grows steeply with max_flips (going from 4 to 8 on a 28-bit/819-barcode codebook was already a
    // ~50x slowdown), which is why misread.qc()/sim.benchmark() in R now default max_flips to 6
    // rather than relying on this fallback (see their own @param max_flips for the reasoning).
    if (max_flips < 1) {max_flips = cb.N_bits;}
    
    // Get info from bc_count_data
    int             N_barcodes    = bc_count_data.nrow();
    int             rate_col      = -1;
    int             variance_col  = -1;
    int             count_col     = -1;
    int             erc_col       = -1; // optional: observed exact-match (uncorrected) read counts
    CharacterVector data_colnames = colnames(bc_count_data);
    for (int i = 0; i < bc_count_data.ncol(); ++i) {
      if (data_colnames[i] == "rates") {rate_col = i;}
      else if (data_colnames[i] == "variance") {variance_col = i;}
      else if (data_colnames[i] == "counts") {count_col = i;}
      else if (data_colnames[i] == "counts_uncorrected") {erc_col = i;}
    }
    if (rate_col == -1 || variance_col == -1 || count_col == -1) {
      Rcpp::stop("bc_count_data must have columns named 'rates', 'variance', and 'counts'.");
    }
   
    // Extract rate, variance, and count data into vectors
    Vdbl bc_rates(N_barcodes);
    Vdbl bc_variance(N_barcodes);
    Vint bc_counts(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      bc_rates[i]    = static_cast<double>(bc_count_data(i, rate_col));
      bc_variance[i] = static_cast<double>(bc_count_data(i, variance_col));
      bc_counts[i]   = static_cast<int>(bc_count_data(i, count_col));
    }
   
    // Build correction table
    LUT correction_table = build_correction_table(cb, maxHam);
    if (verbose) {
      Rcpp::Rcout << "Correction table built with max-correctable Hamming distance " << maxHam << "." << std::endl;
      Rcpp::Rcout << "Correction table size: " << correction_table.size() << " correctable barcodes." << std::endl;
    }
   
    // Invert correction table
    LUTr correction_table_inverted = invert_lookup_table(correction_table);
    if (verbose) {
      Rcpp::Rcout << "Correction table inverted." << std::endl;
    }
   
    // Assemble parsed ST data
    ST_data STdata = {
      bc_rates, bc_variance, bc_counts,
      cb,
      maxHam,
      correction_table, correction_table_inverted,
      EvalResults(), Vint(N_barcodes, 0),           // placeholders
      std::numeric_limits<double>::infinity(),
      1, max_flips, 10                              // n_forks, max_flips, report_freq
    };
   
    // Optional exact-match (uncorrected) read counts for the erc fit term, carried as a column of
    // bc_count_data (rather than a separate argument) so they can never be a different length than
    // the other columns. prepare_STdata_for_mQC() additionally gates has_erc_obs on weights.erc > 0.
    if (erc_col != -1) {
      STdata.has_erc_obs = true;
      STdata.bc_counts_erc_obs.assign(N_barcodes, 0);
      for (int i = 0; i < N_barcodes; ++i) {
        STdata.bc_counts_erc_obs[i] = static_cast<int>(std::round(bc_count_data(i, erc_col)));
      }
    }
   
    return STdata;
   
  }

// Transformation rate: The probability that transform_flips occurs for a spot with true barcode bc, given fr.
double TR(
    const uint64_t   bc,
    const uint64_t   transform_flips,
    const FlipRates& fr,
    double           grad_scale = 1.0,
    Vdbl*            grad       = nullptr,
    Vdbl*            grad2      = nullptr,
    Vdbl*            grad3      = nullptr
  ) {
   
    /*
     * Assumed conditional probability model:
     *   for each bit i (processed in order), the correlation-adjusted flip probability is
     *      p_adj[i] = rate_flip[i] / prod_{j<i, flip[j]=1}(1 - corr(i,j))
     *  The no-flip probability is then 1 - p_adj[i], ensuring the per-bit probabilities sum to 1.
     *  The correction is applied only to the flip branch; the no-flip branch is derived from it.
     *
     * grad (optional): flat parameter gradient vector [rate10 | rate01 | corr1 | corr0]. When non-null,
     *    grad_scale * ∂TR/∂params is accumulated additively into *grad. The caller is responsible for
     *    pre-zeroing *grad before the first accumulating call.
     *
     * grad2 (optional): a second accumulator, added to with the same grad_scale as grad. Used so that
     *    a single backward pass can simultaneously accumulate into both the corrected-count gradient
     *    (always) and the read-count (exact-match) gradient (only for the i where the candidate misread
     *    barcode equals the barcode of interest), without recomputing TR twice.
     *
     * grad3 (optional): a third accumulator, added to with the same grad_scale as grad/grad2. Used by
     *    expected_bc_count() to also accumulate the hit-count gradient (only for the j where the true
     *    barcode equals the barcode of interest) in the same backward pass, again without recomputing TR.
     */
    
    int    N_bits    = fr.rate10.size();
    int    corr_free = N_bits * (N_bits - 1) / 2;
    double log_tr    = 0.0;
    
    // Per-bit quantities: compute log inverse correlation, stored for backward pass.
    Vdbl log_inv_corr(N_bits, 0.0);
    Vdbl log_adj_flip_vec(N_bits);
    Vboo clamped(N_bits, false);
    Vint bit_vec(N_bits), flip_vec(N_bits);
    
    // Forward pass
    for (int i = 0; i < N_bits; ++i) {
      bit_vec[i]  = (bc              >> i) & 1ULL;
      flip_vec[i] = (transform_flips >> i) & 1ULL;
      for (int j = 0; j < i; ++j) {
        if (flip_vec[j]) {
          // Conditional "if (flip_vec[j])" added by Claude Sonnet 4.6; adjustment by log_inv_corr applied only if the flip happens
          int ij = lower_tri_index(i, j);
          log_inv_corr[i] += bit_vec[i] ? std::log(1.0 - fr.corr1[ij]) : std::log(1.0 - fr.corr0[ij]);
        }
      }
      // Log of the correlation-adjusted flip rate for bit i:
      //   log(rate_flip[i]) - sum_{j<i, flip[j]=1} log(1 - corr(i,j))
      double laf = (bit_vec[i] ? fr.log_rate10[i] : fr.log_rate01[i]) - log_inv_corr[i];
      // Guard (added by Claude Sonnet 4.6): clamp adjusted flip probability to (0, 1) so both branches
      //   remain defined. If correlations push the adjusted rate above 1 the model is undefined; 
      //   clamp to the largest representable value strictly below 1 so log1p(-p) stays finite.
      clamped[i]          = (laf >= 0.0);
      log_adj_flip_vec[i] = clamped[i] ? std::log(1.0 - std::numeric_limits<double>::epsilon()) : laf;
      if (flip_vec[i]) {
        log_tr += log_adj_flip_vec[i];
      } else {
        log_tr += std::log1p(-std::exp(log_adj_flip_vec[i]));
      }
    }
    double tr = std::exp(log_tr);
    
    // Backward pass (only when gradient is requested)
    if (grad || grad2 || grad3) {
      for (int i = 0; i < N_bits; ++i) {
        if (clamped[i]) continue;  // δ_i = 0 in the clamped region
        double p_i     = std::exp(log_adj_flip_vec[i]);
        // δ_i = ∂log(TR)/∂log(p_i)
        double delta_i = flip_vec[i] ? 1.0 : (-p_i / (1.0 - p_i));
        double base    = grad_scale * tr * delta_i;
        // Gradient w.r.t. rate_flip[i]: ∂log(p_i)/∂rate_flip = 1/rate_flip
        if (bit_vec[i]) {
          double contrib = base / fr.rate10[i];
          if (grad)  (*grad)[i]  += contrib;
          if (grad2) (*grad2)[i] += contrib;
          if (grad3) (*grad3)[i] += contrib;
        } else {
          double contrib = base / fr.rate01[i];
          if (grad)  (*grad)[N_bits + i]  += contrib;
          if (grad2) (*grad2)[N_bits + i] += contrib;
          if (grad3) (*grad3)[N_bits + i] += contrib;
        }
        // Gradient w.r.t. corr parameters: ∂log(p_i)/∂corr_relevant[k] = 1/(1-corr)
        for (int j = 0; j < i; ++j) {
          if (!flip_vec[j]) continue;  // corr[i,j] only appears when flip_j = 1
          int k = lower_tri_index(i, j);
          if (bit_vec[i]) {
            double contrib = base / (1.0 - fr.corr1[k]);
            if (grad)  (*grad)[2*N_bits + k]  += contrib;
            if (grad2) (*grad2)[2*N_bits + k] += contrib;
            if (grad3) (*grad3)[2*N_bits + k] += contrib;
          } else {
            double contrib = base / (1.0 - fr.corr0[k]);
            if (grad)  (*grad)[2*N_bits + corr_free + k]  += contrib;
            if (grad2) (*grad2)[2*N_bits + corr_free + k] += contrib;
            if (grad3) (*grad3)[2*N_bits + corr_free + k] += contrib;
          }
        }
      }
    }
    
    return tr;
  }

// Function to compute the expected count for a given barcode of interest (BOI)
std::tuple<double, double, double> expected_bc_count(
    const uint64_t   BOI,                     // Barcode of interest (BOI) for which we want to compute expected count after correction
    const FlipRates& fr,                      // FlipRates struct holding rate10, rate01, and corr
    const Vint&      bc_counts,               // Vector of same length as barcodes, giving the ground-truth number of spots with each barcode
    const V64&       barcodes,                // Vector giving all possible true spot barcodes (genes, not blanks)
    const V64&       corrected_to_BOI,        // Vector of barcodes that would be corrected to barcode of interest (BOI)
    int              max_flips,
    Vdbl*            grad_ecc_b = nullptr,    // pre-zeroed flat vector [rate10|rate01|corr1|corr0] of length n_params. When non-null, ∂count_corrected/∂params is accumulated additively into *grad_ecc_b.
    Vdbl*            grad_erc_b = nullptr,    // pre-zeroed flat vector of the same shape as grad_ecc_b. When non-null, ∂count_read/∂params (i.e., the exact-match/uncorrected-read gradient) is accumulated additively into *grad_erc_b.
    Vdbl*            grad_ehc_b = nullptr     // pre-zeroed flat vector of the same shape as grad_ecc_b. When non-null, ∂count_hit/∂params is accumulated additively into *grad_ehc_b.
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
        double scale          = (double)bc_counts[j];
        double tr             = 0.0;
        bool   is_hit_barcode = (barcodes[j] == BOI); // this j's spots are the ones truly of species BOI
        for (int i = 0; i < N_correctable; ++i) {
          // Get bit-flips required to transform the true barcode into the (mis)read barcode
          uint64_t flips_to_BOI = (barcodes[j] ^ corrected_to_BOI[i]) & ((1ULL << N_bits) - 1);
          // Check hamming distance
          if (__builtin_popcountll(flips_to_BOI) > max_flips) {continue;}
          bool   is_exact_match = (corrected_to_BOI[i] == BOI);
          // TR call: also accumulates scale * ∂TR/∂params into grad_ecc_b when non-null, into
          // grad_erc_b when non-null (exact-match term only), and into grad_ehc_b when non-null
          // (true-barcode-is-BOI term only).
          // ... grad_ecc_b tracks ∂count_corrected/∂params because count_corrected = Σ_j Σ_i scale * TR.
          // ... grad_erc_b tracks ∂count_read/∂params because count_read = Σ_j scale * TR (i = exact-match term only).
          // ... grad_ehc_b tracks ∂count_hit/∂params because count_hit = Σ_i scale * TR, for the single j = BOI.
          double tr_ = TR(
            barcodes[j], flips_to_BOI, fr, scale,
            grad_ecc_b,
            is_exact_match ? grad_erc_b : nullptr,
            is_hit_barcode ? grad_ehc_b : nullptr
          );
          tr += tr_;
          if (is_exact_match) {count_read += tr_ * scale;}
        }
        count_corrected += tr * scale;
        if (is_hit_barcode) {count_hit += tr * scale;}
      }
    }
    return {count_read, count_corrected, count_hit};
  }

// Estimate expected barcode counts (read, expected, and hit), as a function of flip rates and true barcode counts, for all barcodes.
std::tuple<Vdbl, Vdbl, Vdbl> expected_bc_counts(
    const FlipRates& fr,                      // FlipRates struct holding rate10, rate01, and corr
    const Vint&      bc_counts,               // Vector of same length as barcodes, giving the ground-truth number of spots with each barcode
    const V64&       barcodes,                // Vector giving all possible true spot barcodes (genes, not blanks)
    const LUTr&      cti,                     // Inverted correction table mapping each barcode index to vector of misread barcodes that would be corrected to it
    int              n_forks,
    int              max_flips,
    Vdbl*            grad_ecc = nullptr,      // pre-zeroed flat matrix stored row-major, size N_barcodes * n_params. When grad_ecc is non-null, ∂ecc[b]/∂params is accumulated into the slice (*grad_ecc)[b*n_params .. (b+1)*n_params-1].
    Vdbl*            grad_erc = nullptr,      // pre-zeroed flat matrix of the same shape as grad_ecc. When non-null (and grad_ecc is also non-null), ∂erc[b]/∂params (the exact-match/uncorrected-read gradient) is accumulated into the analogous slice.
    Vdbl*            grad_ehc = nullptr       // pre-zeroed flat matrix of the same shape as grad_ecc. When non-null (and grad_ecc is also non-null), ∂ehc[b]/∂params (the hit-count gradient, needed for the delta-method PPV CI) is accumulated into the analogous slice. Independent of grad_erc -- either, both, or neither may be requested.
  ) {
    int N_barcodes  = barcodes.size();
    int N_bits      = fr.rate10.size();
    int corr_free   = N_bits * (N_bits - 1) / 2;
    int n_params    = 2*N_bits + 2*corr_free;
    bool do_grad     = (grad_ecc != nullptr);
    bool do_grad_erc = do_grad && (grad_erc != nullptr);
    bool do_grad_ehc = do_grad && (grad_ehc != nullptr);
    Vdbl ecc(N_barcodes, 0.0); // expected corrected counts
    Vdbl erc(N_barcodes, 0.0); // expected read counts
    Vdbl ehc(N_barcodes, 0.0); // expected hit (i.e., correctly read) counts
   
    // Construct batches of barcodes to process in parallel
    std::vector<Vint> barcode_batches(n_forks);
    for (int batch = 0; batch < n_forks; ++batch) {
      int start = (batch * N_barcodes) / n_forks;
      int end   = ((batch + 1) * N_barcodes) / n_forks;
      for (int i = start; i < end; ++i) {barcode_batches[batch].push_back(i);}
    }
    if (n_forks == 1 && barcode_batches[0].size() != N_barcodes) {
      Rcpp::stop("Error in batching barcodes for simulation. Expected all barcodes to be in one batch when n_forks=1.");
    }
   
    // Child payload layout (per barcode in batch):
    //   [erc_0..erc_{N-1} | ecc_0..ecc_{N-1} | ehc_0..ehc_{N-1} | grad_ecc_b0..b{N-1} | grad_erc_b0..b{N-1} (if requested) | grad_ehc_b0..b{N-1} (if requested)]
    // where each grad_*_b is n_params doubles. The grad_ecc block is present when do_grad=true; the
    // grad_erc and grad_ehc blocks are each independently present when do_grad_erc/do_grad_ehc=true,
    // in that fixed order (erc before ehc) whenever both are requested. Run via run_forked_batches()
    // (fork_utils.h): batches run in forked children when n_forks > 1 (POSIX only), else serially
    // in this process -- see that header for why both cases share this same per-batch logic.
    int n_grad_blocks = do_grad ? (1 + (do_grad_erc ? 1 : 0) + (do_grad_ehc ? 1 : 0)) : 0;
   
    run_forked_batches<double>(
      n_forks,
      // batch_payload_size
      [&](int i) -> size_t {
        return static_cast<size_t>(3 + n_grad_blocks * n_params) * barcode_batches[i].size();
      },
      // compute_batch
      [&](int i) -> Vdbl {
        int  N_barcodes_batch = barcode_batches[i].size();
        int  payload_size     = (3 + n_grad_blocks * n_params) * N_barcodes_batch;
        Vdbl erchc_child(payload_size, 0.0);
        Vdbl grad_ecc_b(do_grad     ? n_params : 0, 0.0); // per-barcode ecc-gradient buffer
        Vdbl grad_erc_b(do_grad_erc ? n_params : 0, 0.0); // per-barcode erc-gradient buffer
        Vdbl grad_ehc_b(do_grad_ehc ? n_params : 0, 0.0); // per-barcode ehc-gradient buffer
       
        // Fixed block offsets within the grad region (ecc block always first when do_grad; erc
        // block next if requested; ehc block after that, right after whichever of ecc/erc precedes it)
        int ecc_grad_start = 3 * N_barcodes_batch;
        int erc_grad_start = ecc_grad_start + N_barcodes_batch * n_params;
        int ehc_grad_start = erc_grad_start + (do_grad_erc ? N_barcodes_batch * n_params : 0);
        
        for (int b = 0; b < N_barcodes_batch; ++b) {
          if (do_grad)     std::fill(grad_ecc_b.begin(), grad_ecc_b.end(), 0.0);
          if (do_grad_erc) std::fill(grad_erc_b.begin(), grad_erc_b.end(), 0.0);
          if (do_grad_ehc) std::fill(grad_ehc_b.begin(), grad_ehc_b.end(), 0.0);
          std::tuple<double, double, double> erchc = expected_bc_count(
            barcodes[barcode_batches[i][b]],
            fr, bc_counts, barcodes,
            cti.at(barcode_batches[i][b]),
            max_flips,
            do_grad     ? &grad_ecc_b : nullptr,
            do_grad_erc ? &grad_erc_b : nullptr,
            do_grad_ehc ? &grad_ehc_b : nullptr
          );
          erchc_child[b]                      = std::get<0>(erchc); // expected read count
          erchc_child[b + N_barcodes_batch]   = std::get<1>(erchc); // expected corrected count
          erchc_child[b + 2*N_barcodes_batch] = std::get<2>(erchc); // expected hit count
          if (do_grad) {
            std::copy(grad_ecc_b.begin(), grad_ecc_b.end(),
                      erchc_child.begin() + ecc_grad_start + b*n_params);
          }
          if (do_grad_erc) {
            std::copy(grad_erc_b.begin(), grad_erc_b.end(),
                      erchc_child.begin() + erc_grad_start + b*n_params);
          }
          if (do_grad_ehc) {
            std::copy(grad_ehc_b.begin(), grad_ehc_b.end(),
                      erchc_child.begin() + ehc_grad_start + b*n_params);
          }
        }
        return erchc_child;
      },
      // reduce_batch
      [&](int i, const Vdbl& erchc_child) {
        int N_barcodes_batch = barcode_batches[i].size();

        // Same fixed block offsets as compute_batch used to build this payload (see above)
        int ecc_grad_start = 3 * N_barcodes_batch;
        int erc_grad_start = ecc_grad_start + N_barcodes_batch * n_params;
        int ehc_grad_start = erc_grad_start + (do_grad_erc ? N_barcodes_batch * n_params : 0);

        for (int b = 0; b < N_barcodes_batch; ++b) {
          int bc_idx = barcode_batches[i][b];
          erc[bc_idx] += erchc_child[b];
          ecc[bc_idx] += erchc_child[b + N_barcodes_batch];
          ehc[bc_idx] += erchc_child[b + 2*N_barcodes_batch];
          if (do_grad) {
            std::copy(
              erchc_child.begin() + ecc_grad_start + b*n_params,
              erchc_child.begin() + ecc_grad_start + (b+1)*n_params,
              grad_ecc->begin() + bc_idx * n_params
            );
          }
          if (do_grad_erc) {
            std::copy(
              erchc_child.begin() + erc_grad_start + b*n_params,
              erchc_child.begin() + erc_grad_start + (b+1)*n_params,
              grad_erc->begin() + bc_idx * n_params
            );
          }
          if (do_grad_ehc) {
            std::copy(
              erchc_child.begin() + ehc_grad_start + b*n_params,
              erchc_child.begin() + ehc_grad_start + (b+1)*n_params,
              grad_ehc->begin() + bc_idx * n_params
            );
          }
        }
      },
      // batch_error_context
      [&](int i) -> std::string {
        return "\nBatch size: " + std::to_string(barcode_batches[i].size()) + " barcodes, from " +
               std::to_string(barcode_batches[i][0]) + " to " + std::to_string(barcode_batches[i].back()) + ".";
      }
    );
   
    return {erc, ecc, ehc};
  }

// Estimate misread loss scalar 
double estimate_misread_loss(
    const ST_data& STdata,
    int            max_flips,
    const List&    fliprate_priors
  ) {
    // Grab and set counts ... uses observed corrected counts as ground-truth counts
    Vint BCcounts = STdata.bc_counts;
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

// Estimate the misread-loss scaler and apply it to bc_counts_true (zeroing blanks), in place.
static void apply_misread_loss_rescaling(
    ST_data&    STdata,
    const List& fliprate_priors
  ) {
    double misread_loss_scaler = estimate_misread_loss(STdata, STdata.max_flips, fliprate_priors);
    int    N_barcodes          = STdata.cb.barcodes.size();
    for (int i = 0; i < N_barcodes; ++i) {
      STdata.bc_counts_true[i] = static_cast<int>(std::round(static_cast<double>(STdata.bc_counts[i]) * misread_loss_scaler));
    }
    for (int idx : STdata.cb.blanks) {
      STdata.bc_counts_true[idx] = 0;
    }
   
    /*
     * Shared by prepare_STdata_at() (called with an empty fliprate_priors, matching its
     * pre-existing behavior) and prepare_STdata_for_mQC() (called with the caller's real
     * fliprate_priors).
     */
  }

// Shared setup for mQC() and mQC_msle_at()
static ST_data prepare_STdata_for_mQC(
    NumericMatrix  bc_counts,
    IntegerMatrix  codebook,
    int            maxHam,
    int            n_forks,
    int            max_flips,
    bool           verbose,
    const List&    fliprate_priors,
    const Weights& weights,
    int            max_fs_bc = 40
  ) {
   
    auto STdata = load_STdata(bc_counts, codebook, max_flips, maxHam, verbose);
    STdata.n_forks = n_forks;
   
    // Set up MSLE weights and the flip-rate/correlation prior.
    STdata.msle_weights.assign(STdata.cb.barcodes.size(), 1.0);
    for (int idx : STdata.cb.blanks) {
      STdata.msle_weights[idx] = weights.blank;
    }
    STdata.fr_priors               = pack_fr_priors(fliprate_priors);
    STdata.prior_weight            = weights.prior;
    STdata.dispersion_weight       = weights.dispersion;
    STdata.feasibility_weight      = weights.feasibility;
    STdata.max_fs_bc = max_fs_bc;
   
    // Set up the optional joint erc-fit term: load_STdata() already parsed bc_counts' optional
    // "counts_uncorrected" column into has_erc_obs/bc_counts_erc_obs, if present; just gate it on
    // the fit weight here.
    STdata.erc_weight  = weights.erc;
    STdata.has_erc_obs = STdata.has_erc_obs && weights.erc > 0.0;
   
    // Estimate misread loss scaler and rescale observed counts
    apply_misread_loss_rescaling(STdata, fliprate_priors);
   
    return STdata;
   
    /*
     * Function description:
     *
     * Shared setup for mQC() and mQC_msle_at(): load ST_data and configure everything the penalized MSLE
     * objective (mQC_msle(), via its ST_data* argument) depends on -- msle_weights/fr_priors/prior_weight/
     * dispersion_weight, the optional joint erc-fit term, and the misread-loss rescaling of bc_counts_true.
     * The erc term's observed counts come from bc_counts' optional "counts_uncorrected" column (parsed
     * by load_STdata() into has_erc_obs/bc_counts_erc_obs); here that's just gated on weights.erc > 0, so
     * it can no longer go out of sync with bc_counts the way a separate obs_erc argument could. Callers
     * set n_forks/report_freq themselves (report_freq isn't shared, and n_forks must be set before this
     * runs since estimate_misread_loss() uses it).
     *
     * Initial draft of this function written by Claude Sonnet 5.
     */
  }

/*
 * *********************************************************************************************************************
 * Functions to estimate flip rates, PPV, and CR
 */

// Gaussian-prior + dispersion penalty on flip rates and bit-flip correlations.
double fr_prior_penalty(
    const FlipRates&      fr,
    const FlipRatePriors& priors,
    double                prior_weight,
    double                dispersion_weight,
    bool                  do_grad,
    Vdbl*                 grad = nullptr
  ) {
    int    N_bits    = fr.rate10.size();
    int    corr_free = N_bits * (N_bits - 1) / 2;
    int    n_params  = 2*N_bits + 2*corr_free;
    double total     = 0.0;
   
    // Tiered correlation prior: priors.sdcorr_pairs, if supplied (see build_corr_priors() in R),
    // gives each of the corr_free pairs its own sd instead of sharing the single scalar sdcorr.
    bool has_pair_sd = !priors.sdcorr_pairs.empty();
    if (has_pair_sd && (int)priors.sdcorr_pairs.size() != corr_free) {
      Rcpp::stop("fliprate_priors$sdcorr_pairs must have length N_bits*(N_bits-1)/2 (one sd per bit pair).");
    }
   
    if (prior_weight > 0.0) {
      double penalty = 0.0;
      for (int i = 0; i < N_bits; ++i) {
        double z = (fr.rate10[i] - priors.expected10) / priors.sd10;
        penalty += z * z;
        if (do_grad) (*grad)[i] += prior_weight * 2.0 * z / priors.sd10 / static_cast<double>(n_params);
      }
      for (int i = 0; i < N_bits; ++i) {
        double z = (fr.rate01[i] - priors.expected01) / priors.sd01;
        penalty += z * z;
        if (do_grad) (*grad)[N_bits + i] += prior_weight * 2.0 * z / priors.sd01 / static_cast<double>(n_params);
      }
      for (int i = 0; i < corr_free; ++i) {
        double sdcorr_i = has_pair_sd ? priors.sdcorr_pairs[i] : priors.sdcorr;
        double z1 = (fr.corr1[i] - priors.expectedcorr) / sdcorr_i;
        penalty  += z1 * z1;
        if (do_grad) (*grad)[2*N_bits + i] += prior_weight * 2.0 * z1 / sdcorr_i / static_cast<double>(n_params);
        double z0 = (fr.corr0[i] - priors.expectedcorr) / sdcorr_i;
        penalty  += z0 * z0;
        if (do_grad) (*grad)[2*N_bits + corr_free + i] += prior_weight * 2.0 * z0 / sdcorr_i / static_cast<double>(n_params);
      }
      penalty /= static_cast<double>(n_params);
      total   += prior_weight * penalty;
    }
   
    if (dispersion_weight > 0.0 && corr_free > 1) {
      // When a tiered sdcorr_pairs prior is in use, target the mean of the per-pair variances
      // instead of a single sdcorr^2, so this penalty's "expected heterogeneity" stays consistent
      // with the tiers rather than fighting them (e.g. pulling the near-zero tier back open).
      double target_var;
      if (has_pair_sd) {
        double sumsq = 0.0;
        for (double s : priors.sdcorr_pairs) sumsq += s * s;
        target_var = sumsq / static_cast<double>(priors.sdcorr_pairs.size());
      } else {
        target_var = priors.sdcorr * priors.sdcorr;
      }
     
      auto dispersion_term = [&](const Vdbl& x, int grad_offset) {
        double mean    = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(corr_free);
        double var_emp = 0.0;
        for (double xi : x) var_emp += (xi - mean) * (xi - mean);
        var_emp       /= static_cast<double>(corr_free);
        double diff    = var_emp - target_var;
        if (do_grad) {
          for (int i = 0; i < corr_free; ++i) {
            (*grad)[grad_offset + i] += dispersion_weight * 4.0 * diff * (x[i] - mean) / static_cast<double>(corr_free);
          }
        }
        return diff * diff;
      };
     
      double disp_penalty = dispersion_term(fr.corr1, 2*N_bits) + dispersion_term(fr.corr0, 2*N_bits + corr_free);
      total += dispersion_weight * disp_penalty;
    }
   
    return total;
  }

// nlopt objective callback for L-BFGS minimization of msle.
static double mQC_msle(
    const Vdbl& x,
    Vdbl&       grad,   // nlopt signals "gradient not needed" by passing an empty grad vector
    void*       f_data  // must point to an ST_data struct.
  ) {
    
    // Extract parameters
    ST_data* d          = static_cast<ST_data*>(f_data);
    int      N_bits     = d->cb.N_bits;
    int      N_barcodes = d->cb.barcodes.size();
    int      corr_free  = N_bits * (N_bits - 1) / 2;
    int      n_params   = 2*N_bits + 2*corr_free;
    bool     do_grad    = !grad.empty();
    
    // Sanity check 
    if (x.size() != n_params) {Rcpp::stop("x length doesn't match n_params");}
    
    // Pack flip rates from x
    FlipRates fr = pack_fr(x, N_bits);
    
    // Need gradients? 
    // ... Also need erc's gradient when computing the value-only (do_grad=false) call is fine without it,
    //     but nlopt with LD_LBFGS always requests gradients, so do_erc effectively tracks do_grad here.
    bool do_erc        = do_grad && d->has_erc_obs && d->erc_weight > 0.0;
    bool need_erc_term = d->has_erc_obs && d->erc_weight > 0.0;
    Vdbl grad_ecc_flat, grad_erc_flat;
    if (do_grad) grad_ecc_flat.assign(N_barcodes * n_params, 0.0);
    if (do_erc)  grad_erc_flat.assign(N_barcodes * n_params, 0.0);
    
    // Compute expected counts from flip rates x
    auto erchc = expected_bc_counts(
      fr,
      d->bc_counts_true,
      d->cb.barcodes,
      d->correction_table_inverted,
      d->n_forks,
      d->max_flips,
      do_grad ? &grad_ecc_flat : nullptr,
      do_erc  ? &grad_erc_flat : nullptr
    );
    
    // Compute mean squared log error and weight sum
    double msle       = compute_msle(d->bc_counts, std::get<1>(erchc), d->msle_weights);
    double weight_sum = std::accumulate(d->msle_weights.begin(), d->msle_weights.end(), 0.0);
    
    // Compute gradient of MSLE wrt x
    if (do_grad) {
      std::fill(grad.begin(), grad.end(), 0.0);
      const auto& ecc = std::get<1>(erchc);
      for (int b = 0; b < N_barcodes; ++b) {
        // Chain rule: ∂MSLE/∂θ_k = (w_b/ΣW) * 2*(log(ecc_b+1)-log(obs_b+1))/(ecc_b+1) * ∂ecc_b/∂θ_k
        // ... log equivalent: log(A) - log(B) = log(A/B)
        double coeff = d->msle_weights[b] * 2.0 * std::log((ecc[b] + 1.0) / (static_cast<double>(d->bc_counts[b]) + 1.0)) / (ecc[b] + 1.0) / weight_sum;
        for (int k = 0; k < n_params; ++k) {
          grad[k] += coeff * grad_ecc_flat[b * n_params + k];
        }
      }
    }
    
    // Optional joint loss (and gradient) on exact-match (uncorrected) read counts (erc)
    if (need_erc_term) {
      const auto& erc = std::get<0>(erchc);
      double msle_erc = compute_msle(d->bc_counts_erc_obs, erc, d->msle_weights);
      msle += d->erc_weight * msle_erc;
      if (do_grad) {
        for (int b = 0; b < N_barcodes; ++b) {
          double coeff = d->erc_weight * d->msle_weights[b] * 2.0 * std::log((erc[b] + 1.0) / (static_cast<double>(d->bc_counts_erc_obs[b]) + 1.0)) / (erc[b] + 1.0) / weight_sum;
          for (int k = 0; k < n_params; ++k) {
            grad[k] += coeff * grad_erc_flat[b * n_params + k];
          }
        }
      }
    }
    
    // Gaussian-prior + dispersion penalty on flip rates and bit-flip correlations (shared with the
    // DG fit's MCMCSA objective -- see fr_prior_penalty() above).
    if (d->prior_weight > 0.0 || d->dispersion_weight > 0.0) {
      msle += fr_prior_penalty(fr, d->fr_priors, d->prior_weight, d->dispersion_weight, do_grad, &grad);
    }
   
    // Feasibility penalty: how far corr1/corr0 would have to be distorted to become a valid (PSD)
    // per-barcode Gaussian correlation structure -- see dg_feasibility_penalty() in the "DG
    // simulation functions" section. Keeps the analytic fit's corr1/corr0 comparable to what
    // mQC_DG() could actually simulate, instead of letting them wander somewhere DG would have to
    // heavily correct.
    if (d->feasibility_weight > 0.0) {
      Vdbl feas_grad;
      if (do_grad) {feas_grad.assign(n_params, 0.0);}
      double pen = dg_feasibility_penalty(fr, d->cb, d->max_fs_bc, do_grad, do_grad ? &feas_grad : nullptr);
      msle += d->feasibility_weight * pen;
      if (do_grad) {
        for (int k = 2*N_bits; k < n_params; ++k) {grad[k] += d->feasibility_weight * feas_grad[k];}
      }
    }
   
    // Advance eval counter
    d->eval_results.n_evals++;
    
    // Update eval_results
    if (msle < d->best_msle) {
      int call_n = d->eval_results.n_evals;
      if (call_n == 1 || call_n % d->report_freq == 0) {
        Rcpp::Rcout << "  eval " << call_n << ", msle: " << msle << std::endl;
      }
      // Save eval results 
      d->eval_results.msle_hist.push_back(msle);
      d->eval_results.ehc = std::get<2>(erchc);
      d->eval_results.ecc = std::get<1>(erchc);
      d->eval_results.erc = std::get<0>(erchc);
      // Save best_msle
      d->best_msle = msle;
    }
    
    return msle;
  }

// Main misread quality control function
// [[Rcpp::export]]
List mQC(
    NumericMatrix bc_counts,                   // Optionally include a "counts_uncorrected" column (observed exact-match/uncorrected read counts, same row order as the other columns) to enable the erc loss term; see weights$erc.
    IntegerMatrix codebook,
    int           maxHam,
    int           n_forks,
    int           max_flips,
    int           report_freq     = 1,
    int           maxeval         = 1000,
    List          fliprate_priors = List(),
    List          weights         = List(),   // See Weights/pack_weights(): blank/prior/dispersion/erc/feasibility, all defaulted there
    int           n_restarts      = 1,        // Number of L-BFGS restarts; restart 0 uses the standard init, restarts > 0 randomize the correlation block (1 matches prior behavior)
    int           max_fs_bc       = 40        // Cap on the number of (deterministically subsampled) barcodes used to estimate the feasibility penalty per evaluation, for speed on large codebooks
  ) {
   
    // Load in data and set up the msle objective (weights, priors, erc term, misread-loss rescaling)
    Weights w      = pack_weights(weights);
    auto    STdata = prepare_STdata_for_mQC(
      bc_counts, codebook, maxHam, n_forks, max_flips,
      true, // Make verbose printouts
      fliprate_priors, w,
      max_fs_bc
    );
    STdata.report_freq = report_freq;
   
    // Initialize parameters stochastically from STdata.fr_priors
    int          N_bits        = STdata.cb.N_bits;
    int          N_barcodes    = STdata.cb.barcodes.size(); 
    int          corr_free     = N_bits * (N_bits - 1) / 2;
    size_t       n             = 2*N_bits + 2*corr_free;
    unsigned int mQC_init_seed = 918273645u;
    Vdbl         x0            = draw_random_fr(STdata, STdata.fr_priors, mQC_init_seed);
    Vdbl         x0_           = x0; 
    
    // Set parameter bounds
    double near_one = 1.0 - std::numeric_limits<double>::epsilon(); 
    Vdbl lb(n, std::numeric_limits<double>::epsilon());
    Vdbl ub(n, std::numeric_limits<double>::epsilon()); 
    for (int i = 0; i < N_bits; ++i) {
      ub[i]          =  near_one;
      ub[N_bits + i] =  near_one;
    }
    for (size_t i = 2*N_bits; i < n; ++i) {
      lb[i]          = -near_one;
      ub[i]          =  near_one;
    }
   
    // Initialize matrices to hold results
    NumericMatrix erchc_plus(N_barcodes, 7);
    colnames(erchc_plus) = CharacterVector({"erc", "ecc", "ehc", "CR", "PPV", "PPV_lower", "PPV_upper"});
    NumericVector fr(n);
    
    // Setup nlopt
    double ftol_rel = 1e-8;
    double xtol_rel = 1e-6;
    if (n_restarts < 1) n_restarts = 1;
    boost::random::mt19937 restart_rng(4507919u + static_cast<unsigned int>(N_barcodes) + static_cast<unsigned int>(n));
    double       best_minf = std::numeric_limits<double>::infinity();
    Vdbl         best_x;
    EvalResults  best_eval_results;
    bool         any_success = false;
    Rcpp::Rcout << "\nRunning L-BFGS (nlopt::LD_LBFGS)" << ", maxeval=" << maxeval << ", ftol_rel=" << ftol_rel
                << ", xtol_rel=" << xtol_rel << ", n_restarts=" << n_restarts << std::endl;
    
    // Run L-BFGS via nlopt, with optional multi-start
    for (int r = 0; r < n_restarts; ++r) {
      
      // Set initial parameter (flip rate) values
      Vdbl x0_r = x0;
      if (r > 0) {
        unsigned int seed_r = restart_rng();
        x0_r = draw_random_fr(STdata, STdata.fr_priors, seed_r);
      }
      
      // Reset per-restart tracking so eval_results/best_msle reflect only this restart
      STdata.best_msle    = std::numeric_limits<double>::infinity();
      STdata.eval_results = EvalResults();
      if (n_restarts > 1) Rcpp::Rcout << "\n-- Restart " << (r+1) << "/" << n_restarts << " --" << std::endl;
      
      // Initialize and set optimizer
      double minf_r = std::numeric_limits<double>::quiet_NaN();
      nlopt::opt opt_r(nlopt::LD_LBFGS, n);
      opt_r.set_lower_bounds(lb);
      opt_r.set_upper_bounds(ub);
      opt_r.set_min_objective(mQC_msle, &STdata);
      opt_r.set_ftol_rel(ftol_rel);
      opt_r.set_xtol_rel(xtol_rel);
      opt_r.set_maxeval(maxeval);
      
      // Run optimizer
      try {
        nlopt::result res = opt_r.optimize(x0_r, minf_r);
        Rcpp::Rcout << "  L-BFGS finished (result code " << (int)res << "), final msle: " << minf_r << std::endl;
      } catch (const std::exception& e) {
        Rcpp::Rcout << "  L-BFGS warning: " << e.what() << "\n  Using best value found so far this restart." << std::endl;
        minf_r = STdata.best_msle;
      }
      
      // Check results from this restart and save if best
      if (!std::isnan(minf_r) && minf_r < best_minf) {
        best_minf         = minf_r;
        best_x            = x0_r;
        best_eval_results = STdata.eval_results;
        any_success       = true;
      }
    }
    
    if (!any_success) {
      Rcpp::stop("All restarts failed to produce a finite objective value.");
    }
    if (n_restarts > 1) {
      Rcpp::Rcout << "\nBest of " << n_restarts << " restarts: msle = " << best_minf << std::endl;
    }
    
    x0                  = best_x;
    double minf         = best_minf;
    STdata.eval_results = best_eval_results;
    
    /*
     * ---- Delta-method 95% CI for PPV (Gauss-Newton/Laplace approximation) ----
     * 
     * Var(PPV_b) ~= grad(PPV_b)^T * Sigma_theta * grad(PPV_b), where Sigma_theta is an approximate 
     * parameter covariance estimated from the MSLE term's Gauss-Newton information matrix
     * M = sum_b weight_b * d_b^T d_b, with d_b = d(log(ecc_b+1))/d(theta) = grad_ecc[b,:]/(ecc_b+1)
     * (the same per-barcode chain-rule term mQC_msle()'s own gradient uses). grad(PPV_b) follows from
     * PPV_b = ehc_b/ecc_b via the quotient rule, using the hit-count Jacobian added alongside this CI
     * (see TR()'s grad3 / expected_bc_count()'s grad_ehc_b / expected_bc_counts()'s grad_ehc).
     * 
     * M is built from the MSLE/data term only, deliberately excluding the Gaussian-prior/dispersion
     * penalty curvature: the prior exists to stabilize optimization, not to represent genuine
     * information about a parameter, and folding it in would make the CI look tighter specifically
     * where the *data* is uninformative. This package's benchmarking vignette found most of the
     * bit-flip correlation block to be exactly that (SVD condition number in the billions even at the
     * true parameters) -- M is correspondingly near-singular in those same directions, so its
     * pseudo-inverse is truncated (eigenvalues > 1e-3 * largest kept, matching that same vignette's
     * convention), not inverted outright. The honest consequence is a wide CI in directions the data
     * can't see, rather than a falsely tight one (see the Wald-vs-MCMCSA comparison in
     * old_vs_new_ppv_comparison/, which found exactly this failure mode in the naive single-proportion
     * Wald interval this replaces).
     */
    Vdbl ppv_lower(N_barcodes, NA_REAL), ppv_upper(N_barcodes, NA_REAL);
    {
      int n_params = (int)n;
      Vdbl grad_ecc_flat(N_barcodes * n_params, 0.0);
      Vdbl grad_ehc_flat(N_barcodes * n_params, 0.0);
      FlipRates fr_fit = pack_fr(x0, N_bits);
      expected_bc_counts(
        fr_fit, STdata.bc_counts_true, STdata.cb.barcodes, STdata.correction_table_inverted,
        STdata.n_forks, max_flips, &grad_ecc_flat, nullptr, &grad_ehc_flat
      );
     
      double weight_sum = std::accumulate(STdata.msle_weights.begin(), STdata.msle_weights.end(), 0.0);
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n_params, n_params);
      for (int b = 0; b < N_barcodes; ++b) {
        double ecc_b = STdata.eval_results.ecc[b];
        double w_b   = STdata.msle_weights[b] / weight_sum;
        Eigen::VectorXd d_b(n_params);
        for (int k = 0; k < n_params; ++k) {
          d_b[k] = grad_ecc_flat[b*n_params + k] / (ecc_b + 1.0);
        }
        M.noalias() += w_b * (d_b * d_b.transpose());
      }
     
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
      const Eigen::VectorXd& evals = es.eigenvalues(); // ascending order
      const Eigen::MatrixXd& evecs = es.eigenvectors();
      double lambda_max = evals(n_params - 1);
      double tol        = 1e-3 * lambda_max;
      double sigma2_hat = minf; // final weighted MSLE, used as the common noise-scale estimate
     
      for (int b = 0; b < N_barcodes; ++b) {
        double ecc_b = STdata.eval_results.ecc[b];
        double ehc_b = STdata.eval_results.ehc[b];
        if (ecc_b <= 0.0) continue; // PPV undefined, leave CI as NA
        Eigen::VectorXd grad_ppv(n_params);
        for (int k = 0; k < n_params; ++k) {
          double decc = grad_ecc_flat[b*n_params + k];
          double dehc = grad_ehc_flat[b*n_params + k];
          grad_ppv[k] = (dehc * ecc_b - ehc_b * decc) / (ecc_b * ecc_b);
        }
        double var_ppv = 0.0;
        for (int k = 0; k < n_params; ++k) {
          if (evals(k) <= tol) continue; // truncate near-singular directions
          double proj = evecs.col(k).dot(grad_ppv);
          var_ppv += sigma2_hat * (proj * proj) / evals(k);
        }
        double se  = std::sqrt(std::max(0.0, var_ppv));
        double ppv = ehc_b / ecc_b;
        ppv_lower[b] = std::max(0.0, ppv - 1.96 * se);
        ppv_upper[b] = std::min(1.0, ppv + 1.96 * se);
      }
    }
   
    // Pack results
    for (int k = 0; k < (int)n; ++k) {fr(k)  = x0[k];}
    for (int i = 0; i < N_barcodes; ++i) {
      erchc_plus(i, 0) = STdata.eval_results.erc[i];
      erchc_plus(i, 1) = STdata.eval_results.ecc[i];
      erchc_plus(i, 2) = STdata.eval_results.ehc[i];
      erchc_plus(i, 3) = STdata.eval_results.ecc[i] > 0.0 ? STdata.eval_results.erc[i] / STdata.eval_results.ecc[i] : NA_REAL;
      erchc_plus(i, 4) = STdata.eval_results.ecc[i] > 0.0 ? STdata.eval_results.ehc[i] / STdata.eval_results.ecc[i] : NA_REAL;
      erchc_plus(i, 5) = ppv_lower[i];
      erchc_plus(i, 6) = ppv_upper[i];
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
      _["erctc_plus"]         = erchc_plus, // Values for read count, corrected count, hit count, CR, and PPV analytically implied by fr
      _["msle"]               = minf
    );
  }

/*
 * *********************************************************************************************************************
 * Diagnostic and benchmarking functions related to mQC
 */

// Evaluate mQC()'s penalized MSLE objective at a single, caller-supplied FR, without optimizing.
// [[Rcpp::export]]
double mQC_msle_at(
    const std::vector<double>& FR,
    NumericMatrix              bc_counts,                // Optionally include a "counts_uncorrected" column; see mQC()'s bc_counts.
    IntegerMatrix              codebook,
    int                        maxHam,
    int                        n_forks,
    int                        max_flips,
    List                       fliprate_priors = List(),
    List                       weights         = List(), // See Weights/pack_weights(): blank/prior/dispersion/erc/feasibility, all defaulted there
    int                        max_fs_bc = 40
  ) {

    // Load in data and set up the msle objective, mirroring mQC()'s setup exactly so the result
    // is comparable to mQC()'s "msle"
    Weights w = pack_weights(weights);
    auto STdata = prepare_STdata_for_mQC(
      bc_counts, codebook, maxHam, n_forks, max_flips,
      false, // Suppress verbose printouts
      fliprate_priors, w,
      max_fs_bc
    );
    STdata.report_freq = std::numeric_limits<int>::max();
   
    // Evaluate the penalized MSLE objective at FR (empty grad vector => mQC_msle() skips gradients)
    Vdbl grad;
    Vdbl FR_(FR.begin(), FR.end());
    return mQC_msle(FR_, grad, &STdata);
    
    /*
     * Function description: 
     * 
     * Evaluate the penalized MSLE objective used by mQC() at a single, caller-supplied
     * parameter vector, without running any optimizer. Mirrors mQC()'s STdata setup (msle_weights,
     * fr_priors, prior_weight, dispersion_weight, erc term, misread-loss rescaling of bc_counts_true)
     * exactly, so the returned value is directly comparable to the "msle" mQC() reports for a fitted
     * solution evaluated against the same bc_counts/codebook/weights. Added to let test_fr_recovery()
     * compare the fitted solution's msle against the msle of the true, simulation-stipulated flip
     * rates -- i.e. to tell apart "the optimizer failed to find a better point" from "the penalized
     * objective's optimum isn't at the truth" (see sim.benchmark()'s msle_est/msle_stip output).
     * 
     * Initial draft of this function written by Claude Sonnet 5. 
     */
   
  }

// Shared setup for the single-point diagnostic exports
static std::pair<ST_data, FlipRates> prepare_STdata_at(
    const std::vector<double>& FR,
    NumericMatrix              bc_counts,
    IntegerMatrix              codebook,
    int                        maxHam,
    int                        n_forks,
    int                        max_flips
  ) {
    
    auto STdata = load_STdata(bc_counts, codebook, max_flips, maxHam, false);
    STdata.n_forks = n_forks;
    apply_misread_loss_rescaling(STdata, List());
    FlipRates fr = pack_fr(Vdbl(FR.begin(), FR.end()), STdata.cb.N_bits);
    return {STdata, fr};
   
    /*
     * Function description:
     *
     * Shared setup for the single-point diagnostic exports (mQC_ecc_jacobian_at, mQC_ehc_jacobian_at,
     * mQC_counts_at): load ST_data, estimate and apply the misread-loss rescaling of bc_counts_true exactly
     * as mQC() does (zeroing blanks), then pack FR into a FlipRates struct. These diagnostics skip
     * mQC()'s msle_weights/fr_priors/prior_weight/dispersion_weight/erc setup entirely, since
     * expected_bc_counts() (and hence their gradients/counts) doesn't depend on those -- only
     * mQC_msle()'s penalty terms do.
     * 
     * Initial draft of this function written by Claude Sonnet 5.
     */
  }


// Analytic Jacobian d(ecc_b)/d(theta_k) at a single, caller-supplied FR.
// [[Rcpp::export]]
NumericMatrix mQC_ecc_jacobian_at(
    const std::vector<double>& FR,
    NumericMatrix              bc_counts,
    IntegerMatrix              codebook,
    int                        maxHam,
    int                        n_forks,
    int                        max_flips
  ) {
    
    auto       prepped = prepare_STdata_at(FR, bc_counts, codebook, maxHam, n_forks, max_flips);
    ST_data&   STdata  = prepped.first;
    FlipRates& fr      = prepped.second;
   
    // Compute expected counts with the ecc gradient requested
    int  N_barcodes = STdata.cb.barcodes.size();
    int  n_params   = (int)FR.size();
    Vdbl grad_ecc_flat(N_barcodes * n_params, 0.0);
    
    expected_bc_counts(
      fr,
      STdata.bc_counts_true,
      STdata.cb.barcodes,
      STdata.correction_table_inverted,
      STdata.n_forks,
      max_flips,
      &grad_ecc_flat,
      nullptr
    );
    
    // Unpack the flat row-major (barcode, param) gradient into an N_barcodes x n_params matrix
    NumericMatrix J(N_barcodes, n_params);
    for (int b = 0; b < N_barcodes; ++b) {
      for (int k = 0; k < n_params; ++k) {
        J(b, k) = grad_ecc_flat[b * n_params + k];
      }
    }
    return J;
   
    /*
     * Function description: 
     * 
     * Expose the analytic Jacobian d(ecc_b)/d(theta_k) -- the same gradient expected_bc_counts() already
     * computes internally for L-BFGS, here returned directly instead of being collapsed into a single
     * MSLE gradient. Rows are barcodes, columns are parameters in FR's order (rate10, rate01, corr1, corr0).
     * Added to let identifiability of the bit-flip correlation block be checked directly (SVD of the
     * correlation columns) rather than inferred indirectly from optimizer behavior -- see the cor_recovery
     * diagnostic from sim.benchmark()/analyze_bm(), which found the fitted correlations uncorrelated with
     * the true ones. Mirrors mQC()'s STdata/bc_counts_true setup (misread-loss rescaling) but skips the
     * msle_weights/fr_priors/prior_weight/dispersion_weight setup entirely, since expected_bc_counts() (and
     * hence this Jacobian) doesn't depend on those -- they only enter mQC_msle()'s penalty terms.
     * 
     * Initial draft of this function written by Claude Sonnet 5. 
     */
  }

// Raw (erc, ecc, ehc) count vectors at a single, caller-supplied FR, no gradient.
// [[Rcpp::export]]
NumericMatrix mQC_counts_at(
    const std::vector<double>& FR,
    NumericMatrix              bc_counts,
    IntegerMatrix              codebook,
    int                        maxHam,
    int                        n_forks,
    int                        max_flips
  ) {
   
    auto       prepped = prepare_STdata_at(FR, bc_counts, codebook, maxHam, n_forks, max_flips);
    ST_data&   STdata  = prepped.first;
    FlipRates& fr      = prepped.second;
    int  N_barcodes    = STdata.cb.barcodes.size();
   
    auto erchc = expected_bc_counts(
      fr, STdata.bc_counts_true, STdata.cb.barcodes, STdata.correction_table_inverted,
      STdata.n_forks, max_flips
    );
   
    NumericMatrix out(N_barcodes, 3);
    colnames(out) = CharacterVector({"erc", "ecc", "ehc"});
    for (int b = 0; b < N_barcodes; ++b) {
      out(b, 0) = std::get<0>(erchc)[b];
      out(b, 1) = std::get<1>(erchc)[b];
      out(b, 2) = std::get<2>(erchc)[b];
    }
    return out;
    
    /*
     * Function description: 
     * 
     * Return the raw (erc, ecc, ehc) count vectors at a given FR, no gradient. Added alongside the
     * Jacobian helpers above so their analytic gradients can be checked against finite differences
     * (none of the other exported diagnostics expose count values without also running a fit or a full
     * simulation).
     * 
     * Initial draft of this function written by Claude Sonnet 5. 
     */
  }

// Analytic Jacobian d(ehc_b)/d(theta_k) at a single, caller-supplied FR.
// [[Rcpp::export]]
NumericMatrix mQC_ehc_jacobian_at(
    const std::vector<double>& FR,
    NumericMatrix              bc_counts,
    IntegerMatrix              codebook,
    int                        maxHam,
    int                        n_forks,
    int                        max_flips
  ) {
   
    auto       prepped = prepare_STdata_at(FR, bc_counts, codebook, maxHam, n_forks, max_flips);
    ST_data&   STdata  = prepped.first;
    FlipRates& fr      = prepped.second;
    int  N_barcodes    = STdata.cb.barcodes.size();
    int  n_params      = (int)FR.size();
    Vdbl grad_ecc_flat(N_barcodes * n_params, 0.0); // required alongside grad_ehc (do_grad gates do_grad_ehc)
    Vdbl grad_ehc_flat(N_barcodes * n_params, 0.0);
   
    expected_bc_counts(
      fr,
      STdata.bc_counts_true,
      STdata.cb.barcodes,
      STdata.correction_table_inverted,
      STdata.n_forks,
      max_flips,
      &grad_ecc_flat,
      nullptr,
      &grad_ehc_flat
    );
   
    NumericMatrix J(N_barcodes, n_params);
    for (int b = 0; b < N_barcodes; ++b) {
      for (int k = 0; k < n_params; ++k) {
        J(b, k) = grad_ehc_flat[b * n_params + k];
      }
    }
    return J;
    
    /*
     * Function description: 
     * 
     * Same as mQC_ecc_jacobian_at(), but for the hit-count Jacobian d(ehc_b)/d(theta_k) instead of the
     * corrected-count one. Added for the delta-method PPV confidence interval (PPV_b = ehc_b/ecc_b, so
     * its gradient needs both Jacobians via the quotient rule) -- see mQC()'s CI computation.
     * 
     * Initial draft of this function written by Claude Sonnet 5. 
     */
  }

/*
 * *********************************************************************************************************************
 * Simulation functions
 */

Vint simulate_barcode_spots(
    int                     b,                     // ID of barcode to simulate
    int                     count,                 // Number of spots to simulate
    const FlipRates&        fr,
    const LUT&              correction_table, 
    const V64&              barcodes, 
    boost::random::mt19937& rng
  ) {
    
    // Barcode info
    const int N_barcodes = barcodes.size();
    const int N_bits     = fr.rate10.size(); 
    
    // Vector to hold read, corrected, and hit counts for each barcode
    Vint      rch_counts(3 * N_barcodes, 0); 
    const int read_offset      = 0;
    const int corrected_offset = N_barcodes;
    const int hit_offset       = 2 * N_barcodes;
    
    // Simulate the spots for this barcode, based on TR function
    for (int k = 0; k < count; ++k) {
      
      // Initialize vector to track flips
      Vint flips(N_bits, 0);
      for (int i = 0; i < N_bits; ++i) {
        // Extract bit i of barcode b
        int bit = (barcodes[b] >> i) & 1ULL;
        // Compute correlation adjustment based on previous flips
        double log_inv_corr = 0.0; 
        for (int j = 0; j < i; ++j) {
          if (flips[j]) {
            int ij = lower_tri_index(i, j);
            log_inv_corr += bit ? std::log(1.0 - fr.corr1[ij]) : std::log(1.0 - fr.corr0[ij]);
          }
        }
        // Compute and clamp log of flip rate for this bit
        double FR = (bit ? fr.log_rate10[i] : fr.log_rate01[i]) - log_inv_corr;
        if (FR >= 0.0) {FR = std::log(1.0 - std::numeric_limits<double>::epsilon());}
        // Randomly draw flip or no-flip based on this flip rate
        boost::random::bernoulli_distribution<double> bern(std::exp(FR));
        flips[i] = static_cast<int>(bern(rng));
      }
      // Flip drawn bits in the true barcode
      uint64_t spot_bc_read = barcodes[b] ^ pack(flips);
      
      // Correct decoded label
      auto it              = correction_table.find(spot_bc_read);
      int  label_corrected = it == correction_table.end() ? -1 : it->second;
      int  label_read      = -1;
      
      // ... and correct decoded barcode
      uint64_t spot_bc_corrected = spot_bc_read;
      if (label_corrected >= 0) {
        spot_bc_corrected = barcodes[label_corrected];
        if (spot_bc_corrected == spot_bc_read) {label_read = label_corrected;}
      }
      
      // Add (accumulate) labels
      if (label_read >= 0) {++rch_counts[read_offset + label_read];}
      if (label_corrected >= 0) {
        ++rch_counts[corrected_offset + label_corrected];
        if (label_corrected == b) {++rch_counts[hit_offset + label_corrected];}
      }
      
    }
    
    return rch_counts;
  }

SpotSim make_SpotSim(
    const Vint&      bc_counts,             // Ground-truth counts, per barcode
    const Codebook&  cb, 
    const FlipRates& fr, 
    const LUT&       correction_table,
    int              ran_seed,
    int              n_forks
  ) {
    
    // Initialize spot sim
    SpotSim sim;
    
    // Extract hyperparameters 
    int total_spots = std::accumulate(bc_counts.begin(), bc_counts.end(), 0);
    int N_barcodes  = cb.barcodes.size();
    int N_bits      = cb.N_bits;
    
    // Construct batches of barcodes to simulate in parallel
    std::vector<Vint> barcode_batches(n_forks);
    for (int i = 0; i < N_barcodes; ++i) {barcode_batches[i % n_forks].push_back(i);}
    if (n_forks == 1 && barcode_batches[0].size() != N_barcodes) {
      Rcpp::stop("Error in batching barcodes for simulation. Expected all barcodes to be in one batch when n_forks=1.");
    }
    
    // For each barcode, simulate and decode spot counts
    Rcpp::Rcout << "Simulating and decoding spot counts for each barcode." << std::endl; 
    int cache_size = 3 * N_barcodes;
    Vint rch_counts(cache_size, 0);
    // Run via run_forked_batches() (fork_utils.h): batches run in forked children when
    // n_forks > 1 (POSIX only), else serially in this process.
    run_forked_batches<int>(
      n_forks,
      // batch_payload_size
      [&](int /*i*/) -> size_t {return static_cast<size_t>(cache_size);},
      // compute_batch
      [&](int i) -> Vint {
        Vint rch_counts_batch(cache_size, 0);
        for (int b : barcode_batches[i]) {
          boost::random::mt19937 rng(ran_seed + b);
          Vint temp_vec = simulate_barcode_spots(
            b,
            bc_counts[b],
            fr,
            correction_table,
            cb.barcodes,
            rng
          );
          for (int j = 0; j < cache_size; ++j) {rch_counts_batch[j] += temp_vec[j];}
        }
        return rch_counts_batch;
      },
      // reduce_batch
      [&](int /*i*/, const Vint& temp_vec) {
        for (int j = 0; j < cache_size; ++j) {rch_counts[j] += temp_vec[j];}
      }
    );

    // Resize vectors to count spots per barcode
    sim.read_counts.resize(N_barcodes, 0);
    sim.corrected_counts.resize(N_barcodes, 0);
    sim.hit_counts.resize(N_barcodes, 0);
    
    // Parse out accumulations to spot sim
    const int read_offset      = 0;
    const int corrected_offset = N_barcodes;
    const int hit_offset       = 2 * N_barcodes;
    for (int b = 0; b < N_barcodes; b++) {
      sim.read_counts[b]       = rch_counts[read_offset + b];
      sim.corrected_counts[b]  = rch_counts[corrected_offset + b];
      sim.hit_counts[b]        = rch_counts[hit_offset + b];
    }
    
    return sim;
  }

// Simulate from a "true" FR, fit it back out, and report how well the fit recovers the truth.
// [[Rcpp::export]]
List test_fr_recovery(
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int           n_sims,
    int           maxHam,
    int           n_forks,
    int           max_flips,
    int           report_freq,
    int           maxeval,
    List          fliprate_priors      = List(),
    List          stip_fliprate_priors = List(),
    List          weights              = List(),
    int           n_restarts           = 1
  ) {
    Weights w = pack_weights(weights);
   
    // Load in data
    auto STdata = load_STdata(
      bc_counts, 
      codebook, 
      max_flips,
      maxHam,
      false  // Suppress verbose printouts
    );
    STdata.n_forks     = n_forks;
    STdata.report_freq = report_freq;
    
    // Find bc_counts columns
    int             rate_col      = -1;
    int             variance_col  = -1;
    int             count_col     = -1;
    int             erc_col_in    = -1; // pre-existing "counts_uncorrected" column, if bc_counts already has one
    CharacterVector data_colnames = colnames(bc_counts);
    for (int i = 0; i < bc_counts.ncol(); ++i) {
      if (     data_colnames[i] == "rates")              {rate_col     = i;}
      else if (data_colnames[i] == "variance")            {variance_col = i;}
      else if (data_colnames[i] == "counts")              {count_col    = i;}
      else if (data_colnames[i] == "counts_uncorrected")  {erc_col_in   = i;}
    }
    // ... save original count column
    NumericVector bc_counts_original(bc_counts(_, count_col));
   
    // When the erc term is enabled, work off a copy of bc_counts carrying a "counts_uncorrected"
    // column (filled with each replicate's simulated read counts below) -- reusing one bc_counts
    // already has, or appending one otherwise -- so it can be threaded into mQC()/mQC_msle_at() as a
    // bc_counts column rather than a separate obs_erc argument. Built once here, outside the
    // simulation loop, since rate_col/count_col are progressively rewritten in place across
    // replicates and must stay on the same matrix object. When disabled, bc_counts_fit just aliases
    // bc_counts (no column appended).
    int           erc_col       = -1;
    NumericMatrix bc_counts_fit = bc_counts;
    if (w.erc > 0.0 && erc_col_in != -1) {
      erc_col = erc_col_in; // reuse the existing column; bc_counts_fit already aliases bc_counts
    } else if (w.erc > 0.0) {
      int ncol_orig = bc_counts.ncol();
      erc_col       = ncol_orig;
      bc_counts_fit = NumericMatrix(bc_counts.nrow(), ncol_orig + 1);
      for (int i = 0; i < bc_counts.nrow(); ++i) {
        for (int j = 0; j < ncol_orig; ++j) {bc_counts_fit(i, j) = bc_counts(i, j);}
      }
      CharacterVector colnames_fit(ncol_orig + 1);
      for (int j = 0; j < ncol_orig; ++j) {colnames_fit[j] = data_colnames[j];}
      colnames_fit[ncol_orig] = "counts_uncorrected";
      colnames(bc_counts_fit) = colnames_fit;
    }
   
    // Use observed corrected counts as ground-truth counts for simulation
    Vint BCcounts = STdata.bc_counts;
    // ... with blanks set to zero. 
    for (int i : STdata.cb.blanks) {BCcounts[i] = 0;}
    
    // Draw and pack random flip rates 
    int       N_bits     = STdata.cb.N_bits;
    int       N_barcodes = STdata.cb.barcodes.size(); 
    int       corr_free  = N_bits * (N_bits - 1) / 2;
    size_t    n          = 2*N_bits + 2*corr_free;
    // Truth is drawn from stip_fliprate_priors when supplied, else from fliprate_priors itself
    // (old behavior -- see stip_fliprate_priors' declaration above)
    List      truth_priors = stip_fliprate_priors.size() > 0 ? stip_fliprate_priors : fliprate_priors;
    auto      FR         = draw_random_fr(STdata, pack_fr_priors(truth_priors));
    FlipRates fr         = pack_fr(FR, N_bits);
   
    // Compute expected counts (read, corrected, hit) from simulation ground-truth
    auto expected_in_sim = expected_bc_counts(
      fr, 
      BCcounts, 
      STdata.cb.barcodes, 
      STdata.correction_table_inverted,
      n_forks,
      STdata.max_flips
    );
    Vdbl erc_exp = std::get<0>(expected_in_sim);
    Vdbl ecc_exp = std::get<1>(expected_in_sim);
    Vdbl ehc_exp = std::get<2>(expected_in_sim);
    Vdbl CR__exp(N_barcodes, 0.0);
    Vdbl PPV_exp(N_barcodes, 0.0);
    for (int i = 0; i < N_barcodes; ++i) {
      CR__exp[i] = ecc_exp[i] > 0.0 ? erc_exp[i] / ecc_exp[i] : NA_REAL;
      PPV_exp[i] = ecc_exp[i] > 0.0 ? ehc_exp[i] / ecc_exp[i] : NA_REAL;
    }
    
    // Initialize matrices to hold results 
    // ... flip-rate estimates
    NumericMatrix fr__est(n_sims, n);
    // ... simulated counts
    IntegerMatrix orc_sim(n_sims, N_barcodes);
    IntegerMatrix occ_sim(n_sims, N_barcodes);
    IntegerMatrix ohc_sim(n_sims, N_barcodes);
    NumericMatrix CR__sim(n_sims, N_barcodes); 
    NumericMatrix PPV_sim(n_sims, N_barcodes); 
    // ... analytic estimates based of fits of simulation data
    NumericMatrix erc_est(n_sims, N_barcodes); 
    NumericMatrix ecc_est(n_sims, N_barcodes); 
    NumericMatrix ehc_est(n_sims, N_barcodes); 
    NumericMatrix CR__est(n_sims, N_barcodes);
    NumericMatrix PPV_est(n_sims, N_barcodes);
    // ... penalized-MSLE diagnostic: fitted solution's msle vs. the true/stipulated FR's msle,
    //     both evaluated on the same simulated bc_counts under the same weights/priors, to tell
    //     apart an optimizer failure (msle_stip < msle_est) from a biased objective optimum
    //     (msle_stip >= msle_est)
    NumericVector msle_est(n_sims);
    NumericVector msle_stip(n_sims);
   
    for (int s = 0; s < n_sims; ++s) {
      
      // Simulate corrected reads from this data and stipulated set of flip rates 
      int s_plus = s + 1; 
      Rcpp::Rcout << "\n------ Running Pass Number " << s_plus << "/" << n_sims << " ------" << std::endl;
      Rcpp::Rcout << "\nMaking Bernoulli simulation with given data ..." << std::endl; 
      int ran_seed = 123;
      SpotSim sim = make_SpotSim(
        BCcounts, 
        STdata.cb, 
        fr, 
        STdata.correction_table, 
        ran_seed + s, 
        n_forks
      );
      
      // Use simulation results to rewrite bc_counts_fit (which aliases bc_counts when the erc term
      // is disabled, since no extra column was appended above)
      for (int i = 0; i < N_barcodes; ++i) {
        orc_sim(s, i)                = sim.read_counts[i];
        occ_sim(s, i)                = sim.corrected_counts[i];
        ohc_sim(s, i)                = sim.hit_counts[i];
        CR__sim(s, i)                = static_cast<double>(sim.read_counts[i]) / static_cast<double>(sim.corrected_counts[i]);
        PPV_sim(s, i)                = static_cast<double>(sim.hit_counts[i]) / static_cast<double>(sim.corrected_counts[i]);
        bc_counts_fit(i, rate_col)  /= bc_counts_original(i) / static_cast<double>(std::max(1, sim.corrected_counts[i]));
        bc_counts_fit(i, count_col) = static_cast<double>(sim.corrected_counts[i]);
        if (erc_col != -1) {bc_counts_fit(i, erc_col) = static_cast<double>(sim.read_counts[i]);}
      }
     
      // Run mQC to recover flip rates
      Rcpp::Rcout << "\nRunning misread QC with L-BFGS (nlopt) to recover flip rates...\n" << std::endl;
      // Forward the raw weights list as-is: mQC()'s own pack_weights() unpacks blank/prior/
      // dispersion/erc from it directly, and (since this List never carries a "feasibility" key --
      // sim.benchmark() doesn't expose one) falls back to its own default there, same as before.
      List res = mQC(
        bc_counts_fit, codebook,
        maxHam,
        n_forks, max_flips,
        report_freq, maxeval,
        fliprate_priors,
        weights,
        n_restarts
      );
      NumericMatrix erchc_plus = res["erctc_plus"];
      // ... extract flip-rate vector (size n)
      fr__est(s,_) = Rcpp::as<NumericVector>(res["fliprates"]);
      // ... and count vectors (size N_barcodes)
      erc_est(s,_) = NumericVector(erchc_plus(_, 0));
      ecc_est(s,_) = NumericVector(erchc_plus(_, 1));
      ehc_est(s,_) = NumericVector(erchc_plus(_, 2));
      CR__est(s,_) = NumericVector(erchc_plus(_, 3));
      PPV_est(s,_) = NumericVector(erchc_plus(_, 4));
      msle_est[s]  = Rcpp::as<double>(res["msle"]);
     
      // Evaluate the same penalized objective at the true/stipulated FR, against this replicate's
      // own simulated bc_counts (the exact matrix just handed to mQC() above)
      msle_stip[s] = mQC_msle_at(
        FR, bc_counts_fit, codebook,
        maxHam, n_forks, max_flips,
        fliprate_priors,
        weights
      );
     
    }
    
    return List::create(
      // stipulated flip rates and expected counts, given the stipulated flip rates
      _["fr_stip"] = FR,
      _["erc_exp"] = erc_exp,
      _["ecc_exp"] = ecc_exp,
      _["ehc_exp"] = ehc_exp, 
      _["CR__exp"] = CR__exp,
      _["PPV_exp"] = PPV_exp,
      // simulated counts
      _["orc_sim"] = orc_sim,
      _["occ_sim"] = occ_sim, 
      _["ohc_sim"] = ohc_sim, 
      _["CR__sim"] = CR__sim,
      _["PPV_sim"] = PPV_sim, 
      // ... analytic estimates based of fits of simulation data
      _["fr__est"] = fr__est,
      _["erc_est"] = erc_est,
      _["ecc_est"] = ecc_est,
      _["ehc_est"] = ehc_est,
      _["CR__est"] = CR__est,
      _["PPV_est"] = PPV_est,
      // penalized-objective diagnostic (see msle_est/msle_stip declaration above)
      _["msle_est"]  = msle_est,
      _["msle_stip"] = msle_stip
    );
    
    /*
     * Test 1: By comparing _sim to _exp, can estimate the variation in observations (_sim) and how those
     *  observations compare to the expected counts (_exp). Expect observed observations (_sim) to cluster
     *  around the expected value (_exp).
     *
     * Test 2: By comparing _est to _exp, can estimate how well the analytic model recovers the "true" (expected)
     *  values from (simulated) observations.
     *
     * Test 3: By comparing _sim to _est, can estimate how much the analytic model deviates from a naive count
     *  of observations.
     *
     * stip_fliprate_priors: priors used ONLY for drawing the simulated "true" flip rates/correlations
     * (the FR returned as fr_stip). Defaults to an empty list, which falls back to fliprate_priors
     * itself (old behavior: truth and fit share one prior). Pass a different list -- e.g. from
     * build_corr_priors() with a near-zero sd_other -- to simulate a structured truth while fitting
     * with a different (e.g. unstructured, or differently tiered) prior via fliprate_priors. Lets
     * sim.benchmark() test whether a tiered fitting prior actually helps when the truth really has
     * that structure, versus when it doesn't (see the package's benchmarking notes).
     */
   
  }

/*
 * *********************************************************************************************************************
 * DG simulation functions (mQC_DG)
 */

// Draw n i.i.d. samples from N(0, Sigma) given Sigma's lower-triangular Cholesky factor L.
Eigen::MatrixXd dg_mvn_draw_cholesky(
    const Eigen::MatrixXd&  L,
    int                     n,
    boost::random::mt19937& rng
  ) {
    int dim = L.rows();
    boost::random::normal_distribution<double> std_norm(0.0, 1.0);
    Eigen::MatrixXd Z(n, dim);
    for (int r = 0; r < n; ++r) {
      for (int c = 0; c < dim; ++c) {Z(r, c) = std_norm(rng);}
    }
    return Z * L.transpose(); // row s = L * z_s; mean is always 0 here (thresholds carry bit-value info)
   
    /*
     * Sigma = L*L^T. Native, via a caller-supplied boost::random::mt19937 (NOT R's own RNG -- this
     * is called from inside forked child processes in make_SpotSim_DG() below, and R's
     * embedded-interpreter RNG state is not fork-safe to touch there; see make_SpotSim()'s own
     * boost::random::mt19937 use for precedent). Mathematically identical to MASS::mvrnorm (any
     * valid matrix square root produces the same multivariate normal law) -- this just avoids the
     * R callback dg_legacy.cpp's generate_spots_ci() pays per barcode per cell.
     */
  }

// Compute lower-tail Poisson CDF
double dg_ppois_lower(double x, double lambda) {
    if (x < 0.0) {return 0.0;}
    if (lambda <= 0.0) {return 1.0;}
    double x_floor = std::floor(x + 1e-7);
    return boost::math::cdf(boost::math::poisson_distribution<double>(lambda), x_floor);
    /*
     * R::ppois(x, lambda, lower_tail = 1, log_p = 0) semantics, natively: the lower-tail Poisson CDF is 
     * the regularized upper incomplete gamma function gamma_q(floor(x)+1, lambda) (boost::math::cdf()
     * evaluates gamma_q(k+1, lambda) directly without flooring non-integer k -- see
     * boost/math/distributions/poisson.hpp -- so the floor(x + 1e-7) step below, matching R's own
     * nmath/ppois.c, is reproduced explicitly here). Also reproduces R's edge cases: x < 0 -> 0,
     * lambda == 0 -> 1 (boost's cdf() returns 0 for mean == 0, since it has no such Poisson(0)-is-a-
     * point-mass-at-0 special case).
     */
  }

// This step's total "true" transcript counts per barcode, for a freshly-simulated population of N_cells.
Vint dg_draw_true_counts(
    const ST_data&          STdata,
    const FlipRates&        fr,
    int                     N_cells,
    boost::random::mt19937& rng
  ) {
    int N_bits     = STdata.cb.N_bits;
    int N_barcodes = STdata.cb.barcodes.size();

    double mean_Hamming_weight = 0.0;
    for (uint64_t bc : STdata.cb.barcodes) {mean_Hamming_weight += __builtin_popcountll(bc);}
    mean_Hamming_weight /= static_cast<double>(N_barcodes);

    double mean_rate10 = std::accumulate(fr.rate10.begin(), fr.rate10.end(), 0.0) / static_cast<double>(N_bits);
    double mean_rate01 = std::accumulate(fr.rate01.begin(), fr.rate01.end(), 0.0) / static_cast<double>(N_bits);
    double expected_flip_rate =
      (mean_rate10 * mean_Hamming_weight + mean_rate01 * (N_bits - mean_Hamming_weight)) / static_cast<double>(N_bits);
    double expected_decoding_rate = dg_ppois_lower(STdata.maxHam / 2.0, expected_flip_rate * N_bits);
    if (expected_decoding_rate <= 0.0) {expected_decoding_rate = 1e-6;} // guard against pathological flip rates

    Vint true_counts(N_barcodes, 0); // blanks stay 0
    for (int idx : STdata.cb.genes) {
      double gene_kernel_rate = STdata.bc_rates[idx] / expected_decoding_rate;
      double gamma_variance   = STdata.bc_variance[idx];
      double shape = gene_kernel_rate * gene_kernel_rate / gamma_variance;
      double scale = gamma_variance / gene_kernel_rate;
      int    total = 0;
      for (int c = 0; c < N_cells; ++c) {
        double lambda = boost::random::gamma_distribution<double>(shape, scale)(rng);
        total += boost::random::poisson_distribution<int>(lambda)(rng);
      }
      true_counts[idx] = total;
    }
    return true_counts;
   
    /*
     * Via the same over-dispersed-Poisson (Gamma-Poisson) model .DG_set_true_spot_info() in
     * R/dg_legacy.R uses: per cell, draw a Gamma-distributed transcription rate per gene
     * (capturing cell-level over-dispersion), then a Poisson count from that rate. Blanks are
     * always 0 (no true transcripts; their observed counts arise entirely from misreads).
     */
  }

// Simulate `count` spots truly of barcode b and decode/correct them.
Vint simulate_barcode_spots_DG(
    int                     b,
    int                     count,
    const Vdbl&             thresholds,
    const Eigen::MatrixXd&  L,
    const LUT&              correction_table,
    const V64&              barcodes,
    int                     N_bits,
    boost::random::mt19937& rng
  ) {
    int  N_barcodes = barcodes.size();
    Vint rch_counts(3 * N_barcodes, 0);
    if (count <= 0) {return rch_counts;}
    const int read_offset      = 0;
    const int corrected_offset = N_barcodes;
    const int hit_offset       = 2 * N_barcodes;
   
    Eigen::MatrixXd draws = dg_mvn_draw_cholesky(L, count, rng);
   
    for (int s = 0; s < count; ++s) {
      Vint read_bits(N_bits);
      for (int i = 0; i < N_bits; ++i) {read_bits[i] = draws(s, i) > thresholds[i] ? 1 : 0;}
      uint64_t spot_bc_read = pack(read_bits);
     
      auto it              = correction_table.find(spot_bc_read);
      int  label_corrected = it == correction_table.end() ? -1 : it->second;
      int  label_read      = -1;
      if (label_corrected >= 0) {
        uint64_t spot_bc_corrected = barcodes[label_corrected];
        if (spot_bc_corrected == spot_bc_read) {label_read = label_corrected;}
      }
      if (label_read >= 0) {++rch_counts[read_offset + label_read];}
      if (label_corrected >= 0) {
        ++rch_counts[corrected_offset + label_corrected];
        if (label_corrected == b) {++rch_counts[hit_offset + label_corrected];}
      }
    }
    return rch_counts;
   
    /*
     * Given this barcode's precomputed per-bit thresholds and Gaussian-correlation Cholesky
     * factor (both built once, by dg_precompute_barcode_noise() below). Returns the same
     * [read|corrected|hit] counts-per-barcode layout make_SpotSim()/simulate_barcode_spots() use.
     */
  }

// Build barcode b's per-bit thresholds and Gaussian-correlation Cholesky factor.
std::pair<Vdbl, Eigen::MatrixXd> dg_precompute_barcode_noise(
    int              b,
    const FlipRates& fr,
    const V64&       barcodes,
    int              N_bits
  ) {
    Vint barcode_bits(N_bits);
    Vdbl thresholds(N_bits);
    for (int i = 0; i < N_bits; ++i) {
      barcode_bits[i] = (barcodes[b] >> i) & 1ULL;
      thresholds[i]   = dg_bit_threshold(barcode_bits[i], fr.rate10[i], fr.rate01[i]);
    }
    Eigen::MatrixXd raw       = dg_raw_gaussian_corr(barcode_bits, fr, thresholds);
    Eigen::MatrixXd projected = dg_nearest_corr_matrix(raw);
    Eigen::LLT<Eigen::MatrixXd> llt(projected);
    Eigen::MatrixXd L = (llt.info() == Eigen::Success)
      ? Eigen::MatrixXd(llt.matrixL())
      : Eigen::MatrixXd::Identity(N_bits, N_bits); // pathological fallback: treat as uncorrelated
    return {thresholds, L};
   
    /*
     * Raw matrix -> nearPD projection -> Cholesky. Pure C++ (dg_nearest_corr_matrix() is now a
     * native nearPD port, no R callback -- see dg_nearPD_native()); split out from
     * simulate_barcode_spots_DG() per barcode so the per-call noise model is built once per active
     * barcode rather than once per simulated spot.
     */
  }

// Simulate and decode spot counts for every barcode with a nonzero true count this step, via DG spot generation.
SpotSim make_SpotSim_DG(
    const Vint&      true_counts,
    const Codebook&  cb,
    const FlipRates& fr,
    const LUT&       correction_table,
    int              ran_seed,
    int              n_forks
  ) {
    SpotSim sim;
    int N_barcodes = cb.barcodes.size();
    int N_bits     = cb.N_bits;
   
    // Active barcodes (nonzero true count this step). Each active barcode's noise model (nearPD
    // projection + Cholesky factor, via dg_precompute_barcode_noise()) is computed inside whichever
    // process -- a forked child, or this process itself when n_use==1 -- ends up simulating that
    // barcode's spots, rather than upfront here for all active barcodes; see this function's closing
    // comment for why that's both correct and fork-safe.
    Vint active;
    for (int b = 0; b < N_barcodes; ++b) {if (true_counts[b] > 0) {active.push_back(b);}}

    // Batch active barcodes across forks (same pattern as make_SpotSim())
    int n_active = active.size();
    int n_use    = std::max(1, std::min(n_forks, std::max(1, n_active)));
    std::vector<Vint> barcode_batches(n_use);
    for (int i = 0; i < n_active; ++i) {barcode_batches[i % n_use].push_back(active[i]);}
   
    int  cache_size = 3 * N_barcodes;
    Vint rch_counts(cache_size, 0);
   
    // Run via run_forked_batches() (fork_utils.h): batches run in forked children when
    // n_use > 1 (POSIX only), else serially in this process.
    run_forked_batches<int>(
      n_use,
      // batch_payload_size
      [&](int /*i*/) -> size_t {return static_cast<size_t>(cache_size);},
      // compute_batch
      [&](int i) -> Vint {
        Vint rch_counts_batch(cache_size, 0);
        for (int b : barcode_batches[i]) {
          auto                    pre = dg_precompute_barcode_noise(b, fr, cb.barcodes, N_bits);
          boost::random::mt19937  rng(ran_seed + b);
          Vint temp_vec = simulate_barcode_spots_DG(
            b, true_counts[b], pre.first, pre.second, correction_table, cb.barcodes, N_bits, rng
          );
          for (int j = 0; j < cache_size; ++j) {rch_counts_batch[j] += temp_vec[j];}
        }
        return rch_counts_batch;
      },
      // reduce_batch
      [&](int /*i*/, const Vint& temp_vec) {
        for (int j = 0; j < cache_size; ++j) {rch_counts[j] += temp_vec[j];}
      }
    );

    sim.read_counts.resize(N_barcodes, 0);
    sim.corrected_counts.resize(N_barcodes, 0);
    sim.hit_counts.resize(N_barcodes, 0);
    const int read_offset      = 0;
    const int corrected_offset = N_barcodes;
    const int hit_offset       = 2 * N_barcodes;
    for (int b = 0; b < N_barcodes; b++) {
      sim.read_counts[b]      = rch_counts[read_offset + b];
      sim.corrected_counts[b] = rch_counts[corrected_offset + b];
      sim.hit_counts[b]       = rch_counts[hit_offset + b];
    }
    return sim;
   
    /*
     * Via DG (correlated-Gaussian-noise) spot generation. Mirrors make_SpotSim()'s fork/pipe
     * parallelism (batches barcodes across n_forks processes, same pipe-based int-buffer handoff).
     * Each active barcode's precompute (dg_precompute_barcode_noise(): nearPD projection + Cholesky
     * -- the heavier numerical cost per active barcode, an iterative eigendecomposition-based
     * projection) now runs inside whichever process actually simulates that barcode, rather than
     * serially in the parent for every active barcode before forking. That used to be forced by
     * dg_nearest_corr_matrix()'s R callback (not fork-safe); now that it's a native nearPD port
     * (dg_nearPD_native()), it's itself fork-safe, and distributing it across n_forks this way needs
     * no new pipe payload -- each child only ever needs the precompute results for barcodes it alone
     * simulates, and (as before) sends back just its simulated counts.
     */
  }

// Median of x, skipping NaN entries (NA_REAL if nothing finite is left), matching R's median().
double dg_median_skipna(Vdbl x) {
  Vdbl finite;
  finite.reserve(x.size());
  for (double xi : x) {if (!std::isnan(xi)) {finite.push_back(xi);}}
  if (finite.empty()) {return NA_REAL;}
  std::sort(finite.begin(), finite.end());
  size_t n = finite.size();
  return n % 2 == 1 ? finite[n / 2] : (finite[n / 2 - 1] + finite[n / 2]) / 2.0;

  /* Matches R's median() with NAs removed: for an even count of finite values, the mean of the
   * two middle ones. */
}

// Quantile of x at probability q, skipping NaN entries, matching R's default quantile(type = 7).
double dg_quantile_skipna(Vdbl x, double q) {
  Vdbl finite;
  finite.reserve(x.size());
  for (double xi : x) {if (!std::isnan(xi)) {finite.push_back(xi);}}
  if (finite.empty()) {return NA_REAL;}
  std::sort(finite.begin(), finite.end());
  size_t n = finite.size();
  if (n == 1) {return finite[0];}
  double h       = (static_cast<double>(n) - 1.0) * q;
  size_t lo      = static_cast<size_t>(std::floor(h));
  size_t hi      = static_cast<size_t>(std::ceil(h));
  return finite[lo] + (h - static_cast<double>(lo)) * (finite[hi] - finite[lo]);
}

DGStepResult dg_simulation_step(
    const ST_data&          STdata,
    const FlipRates&        fr,
    int                     N_cells,
    double                  blank_weight,
    bool                    maintain_gene_identity,
    int                     ran_seed,
    int                     n_forks,
    boost::random::mt19937& rng
  ) {
    const Codebook& cb       = STdata.cb;
    int             N_genes  = cb.genes.size();
    int             N_blanks = cb.blanks.size();

    // 1. Draw this step's true counts and simulate/decode spots
    Vint true_counts = dg_draw_true_counts(STdata, fr, N_cells, rng);
    SpotSim sim = make_SpotSim_DG(true_counts, cb, fr, STdata.correction_table, ran_seed, n_forks);
   
    // 2. PPV per gene (hit / corrected; NaN when corrected == 0)
    Vdbl PPV_genes(N_genes, std::numeric_limits<double>::quiet_NaN());
    for (int g = 0; g < N_genes; ++g) {
      int idx = cb.genes[g];
      if (sim.corrected_counts[idx] > 0) {
        PPV_genes[g] = static_cast<double>(sim.hit_counts[idx]) / static_cast<double>(sim.corrected_counts[idx]);
      }
    }
    
    // 3. Sorted-count comparison (observed vs. simulated corrected counts), weighted log-MSE --
    // drives the MCMCSA acceptance criterion (sim_mse); the sorted arrays themselves aren't kept.
    Vdbl obs_counts_genes(N_genes), obs_counts_blanks(N_blanks);
    for (int g = 0; g < N_genes;  ++g) {obs_counts_genes[g]  = STdata.bc_rates[cb.genes[g]]  * N_cells;}
    for (int b = 0; b < N_blanks; ++b) {obs_counts_blanks[b] = STdata.bc_rates[cb.blanks[b]] * N_cells;}
    Vint obs_gene_order(N_genes), obs_blank_order(N_blanks);
    for (int g = 0; g < N_genes;  ++g) {obs_gene_order[g]  = g;}
    for (int b = 0; b < N_blanks; ++b) {obs_blank_order[b] = b;}
    std::sort(obs_gene_order.begin(),  obs_gene_order.end(),  [&](int a, int b){return obs_counts_genes[a]  > obs_counts_genes[b];});
    std::sort(obs_blank_order.begin(), obs_blank_order.end(), [&](int a, int b){return obs_counts_blanks[a] > obs_counts_blanks[b];});
   
    Vdbl sim_counts_genes(N_genes), sim_counts_blanks(N_blanks);
    for (int g = 0; g < N_genes;  ++g) {sim_counts_genes[g]  = static_cast<double>(sim.corrected_counts[cb.genes[g]]);}
    for (int b = 0; b < N_blanks; ++b) {sim_counts_blanks[b] = static_cast<double>(sim.corrected_counts[cb.blanks[b]]);}
    
    Vdbl counts_sorted(N_genes + N_blanks), obs_counts_sorted(N_genes + N_blanks);
    for (int g = 0; g < N_genes;  ++g) {obs_counts_sorted[g]           = obs_counts_genes[obs_gene_order[g]];}
    for (int b = 0; b < N_blanks; ++b) {obs_counts_sorted[N_genes + b] = obs_counts_blanks[obs_blank_order[b]];}
   
    if (maintain_gene_identity) {
      for (int g = 0; g < N_genes;  ++g) {counts_sorted[g]           = sim_counts_genes[obs_gene_order[g]];}
      for (int b = 0; b < N_blanks; ++b) {counts_sorted[N_genes + b] = sim_counts_blanks[obs_blank_order[b]];}
    } else {
      Vint sim_gene_order(N_genes), sim_blank_order(N_blanks);
      for (int g = 0; g < N_genes;  ++g) {sim_gene_order[g]  = g;}
      for (int b = 0; b < N_blanks; ++b) {sim_blank_order[b] = b;}
      std::sort(sim_gene_order.begin(),  sim_gene_order.end(),  [&](int a, int b){return sim_counts_genes[a]  > sim_counts_genes[b];});
      std::sort(sim_blank_order.begin(), sim_blank_order.end(), [&](int a, int b){return sim_counts_blanks[a] > sim_counts_blanks[b];});
      for (int g = 0; g < N_genes;  ++g) {counts_sorted[g]           = sim_counts_genes[sim_gene_order[g]];}
      for (int b = 0; b < N_blanks; ++b) {counts_sorted[N_genes + b] = sim_counts_blanks[sim_blank_order[b]];}
    }
   
    double gene_weight = 1.0 + (static_cast<double>(N_blanks) / static_cast<double>(N_genes)) * (1.0 - blank_weight);
    double sim_mse = 0.0;
    for (int i = 0; i < N_genes + N_blanks; ++i) {
      double w  = i < N_genes ? gene_weight : blank_weight;
      double se = std::pow(std::log(counts_sorted[i] + 1.0) - std::log(obs_counts_sorted[i] + 1.0), 2.0);
      sim_mse  += se * w;
    }
    sim_mse /= static_cast<double>(N_genes + N_blanks);
   
    DGStepResult out;
    out.sim_mse          = sim_mse;
    out.PPV_genes        = PPV_genes;
    out.read_counts      = Vdbl(sim.read_counts.begin(), sim.read_counts.end());
    out.corrected_counts = Vdbl(sim.corrected_counts.begin(), sim.corrected_counts.end());
    out.hit_counts       = Vdbl(sim.hit_counts.begin(), sim.hit_counts.end());
    return out;
  }

// Draw k distinct indices from {0,...,n-1} without replacement.
Vint dg_sample_without_replacement(
    int                     n,
    int                     k,
    boost::random::mt19937& rng
  ) {
    Vint idx(n);
    for (int i = 0; i < n; ++i) {idx[i] = i;}
    int k_use = std::max(0, std::min(k, n));
    for (int i = 0; i < k_use; ++i) {
      int j = boost::random::uniform_int_distribution<int>(i, n - 1)(rng);
      std::swap(idx[i], idx[j]);
    }
    idx.resize(k_use);
    return idx;

    /*
     * Partial Fisher-Yates via the caller's boost::random::mt19937, so the whole MCMCSA walk's
     * reproducibility is controlled by mQC_DG()'s own `seed` argument, not R's set.seed() -- matching
     * dg_draw_true_counts()'s rng threading and this file's general boost::random convention. Only
     * ever called from the parent process (the outer MCMCSA loop in mQC_DG(), never from inside
     * make_SpotSim_DG()'s forked children).
     */
  }

// C++ engine behind misread.qc.DG(): fits FlipRates via Markov Chain Monte Carlo Simulated Annealing.
// [[Rcpp::export]]
List mQC_DG(
    NumericMatrix bc_counts,
    IntegerMatrix codebook,
    int           maxHam,
    int           resamples,
    int           burnin,
    int           FR_search,
    double        max_step_size,
    int           min_step_size_divider,
    double        max_temperature,
    double        min_temperature_scale,
    int           n_tracker_updates,
    int           n_forks,
    int           N_cells,
    List          fliprate_priors          = List(),
    List          weights                  = List(),
    int           max_fs_bc = 40,
    int           max_calls_per_step = 0,
    unsigned int  seed = 726354819
  ) {

    ST_data STdata = load_STdata(bc_counts, codebook, codebook.ncol(), maxHam, true);
    STdata.n_forks = n_forks;
    FlipRatePriors fr_priors = pack_fr_priors(fliprate_priors);
    Weights        w         = pack_weights(weights);

    int N_bits     = STdata.cb.N_bits;
    int N_genes    = STdata.cb.genes.size();
    int N_blanks   = STdata.cb.blanks.size();
    int corr_free  = N_bits * (N_bits - 1) / 2;
    int n_params   = 2*N_bits + 2*corr_free;
    if (N_genes == 0 || N_blanks == 0) {
      Rcpp::stop("Codebook must contain at least one gene and one blank barcode.");
    }

    double min_step_size   = max_step_size / min_step_size_divider;
    double min_temperature = max_temperature / (min_step_size_divider / min_temperature_scale);

    auto lin_seq = [](double from, double to, int n) {
      Vdbl out(n);
      if (n == 1) {out[0] = from; return out;}
      for (int i = 0; i < n; ++i) {out[i] = from + (to - from) * static_cast<double>(i) / static_cast<double>(n - 1);}
      return out;
    };
    auto lin_seq_int = [](int from, int to, int n) {
      Vint out(n);
      if (n == 1) {out[0] = from; return out;}
      for (int i = 0; i < n; ++i) {
        out[i] = static_cast<int>(std::round(from + (to - from) * static_cast<double>(i) / static_cast<double>(n - 1)));
      }
      return out;
    };

    Vdbl step_size_schedule   = lin_seq(max_step_size, min_step_size, FR_search);
    Vdbl temperature_schedule = lin_seq(max_temperature, min_temperature, FR_search);
    Vdbl blank_weight_schedule = lin_seq(w.max_blank, w.min_blank, FR_search);
    Vint n_retained_bits       = lin_seq_int(N_bits - 1, 0, FR_search);
    Vint n_retained_corrs      = lin_seq_int(corr_free - 1, 0, FR_search);
    {
      Vdbl s2 = lin_seq(max_step_size, min_step_size, burnin);
      Vdbl t2 = lin_seq(max_temperature, min_temperature, burnin);
      Vdbl b2 = lin_seq(w.max_blank, w.min_blank, burnin);
      Vint rb2 = lin_seq_int(N_bits - 1, 0, burnin);
      Vint rc2 = lin_seq_int(corr_free - 1, 0, burnin);
      step_size_schedule.insert(step_size_schedule.end(), s2.begin(), s2.end());
      temperature_schedule.insert(temperature_schedule.end(), t2.begin(), t2.end());
      blank_weight_schedule.insert(blank_weight_schedule.end(), b2.begin(), b2.end());
      n_retained_bits.insert(n_retained_bits.end(), rb2.begin(), rb2.end());
      n_retained_corrs.insert(n_retained_corrs.end(), rc2.begin(), rc2.end());
    }
    step_size_schedule.insert(step_size_schedule.end(), resamples, min_step_size);
    temperature_schedule.insert(temperature_schedule.end(), resamples, min_temperature);
    blank_weight_schedule.insert(blank_weight_schedule.end(), resamples, w.min_blank);
    n_retained_bits.insert(n_retained_bits.end(), resamples, 0);
    n_retained_corrs.insert(n_retained_corrs.end(), resamples, 0);

    // Per-pair proposal-step-size tiering: tightly-pinned prior tiers (small sdcorr_pairs) get
    // smaller per-step Gaussian proposals, not just a bigger acceptance penalty.
    Vdbl corr_step_scale(corr_free, 1.0);
    if (!fr_priors.sdcorr_pairs.empty() && fr_priors.sdcorr > 0.0) {
      for (int k = 0; k < corr_free; ++k) {corr_step_scale[k] = fr_priors.sdcorr_pairs[k] / fr_priors.sdcorr;}
    }

    int N_barcodes = STdata.cb.barcodes.size();
    int n_steps = FR_search + burnin + resamples;
    NumericVector msle(n_steps);
    NumericMatrix PPV_genes_mat(n_steps, N_genes);
    Vdbl stip_rate10_hist(n_steps, NA_REAL), stip_rate01_hist(n_steps, NA_REAL);

    // Per-barcode read/corrected/hit counts, kept only for accepted steps in the post-burn-in
    // resamples window and never themselves sent to R -- collapsed into the small erctc_plus matrix
    // (N_barcodes x 7) once the walk finishes, so the only large structure ever built here stays
    // internal to C++.
    int resample_window_start = FR_search + burnin; // 1-indexed step count at which the window begins
    std::vector<Vdbl> resample_read, resample_corrected, resample_hit;
    resample_read.reserve(resamples + 1);
    resample_corrected.reserve(resamples + 1);
    resample_hit.reserve(resamples + 1);

    // Per-bit/per-pair flip-rate and correlation state, same resamples window, collapsed to a single
    // "fliprates" vector (median per parameter) once the walk finishes -- misread.qc.DG()'s analogue
    // of misread.qc()'s own point-estimate fliprates vector, same length/layout (2*N_bits rates then
    // 2*corr_free correlations), just a posterior median across the chain instead of an L-BFGS optimum.
    std::vector<Vdbl> resample_rate10, resample_rate01, resample_corr1, resample_corr0;
    resample_rate10.reserve(resamples + 1);
    resample_rate01.reserve(resamples + 1);
    resample_corr1.reserve(resamples + 1);
    resample_corr0.reserve(resamples + 1);

    boost::random::mt19937 rng(seed);
    boost::random::normal_distribution<double> std_normal(0.0, 1.0);
    auto fresh_seed = [&rng]() {return boost::random::uniform_int_distribution<unsigned int>()(rng);};

    Vdbl FR_current = draw_random_fr(STdata, fr_priors, fresh_seed());
    FlipRates fr_current = pack_fr(FR_current, N_bits);

    Rcpp::Rcout << "\n\nRunning initial simulation with structured-prior-drawn flip rates and correlations" << std::endl;
    DGStepResult init_result = dg_simulation_step(STdata, fr_current, N_cells, w.max_blank, false, (int)fresh_seed(), n_forks, rng);
    double sim_mse_current = init_result.sim_mse;
    if (w.prior > 0.0 || w.dispersion > 0.0) {
      sim_mse_current += fr_prior_penalty(fr_current, fr_priors, w.prior, w.dispersion, false, nullptr);
    }
    if (w.feasibility > 0.0) {
      sim_mse_current += w.feasibility * dg_feasibility_penalty(fr_current, STdata.cb, max_fs_bc);
    }
    Rcpp::Rcout << "\nInitial mse: " << std::round(sim_mse_current * 1000.0) / 1000.0 << std::endl;

    Rcpp::Rcout << "\n\nRunning Markov chain steps: \n\nRunning initial steps ..." << std::endl;
    int    step = 0, ctr = 0, last_step = 0, last_ctr = 0, last_message = 0;
    int    attempts_this_step = 0; // resets whenever step advances; see max_calls_per_step above
    bool   maintain_gene_identity = false;
    double rate10_mean = NA_REAL, rate01_mean = NA_REAL;
    int    FR_sample_start = static_cast<int>(FR_search * 0.8);
    Vint   progress_tracker = Vint();
    {
      Vdbl pt = lin_seq(1.0, static_cast<double>(n_steps), n_tracker_updates + 1);
      for (double p : pt) {progress_tracker.push_back(static_cast<int>(std::round(p)));}
      progress_tracker.back() = n_steps;
    }

    while (step < n_steps) {

      ctr++;
      attempts_this_step++;

      Vdbl rate10_next(N_bits), rate01_next(N_bits);
      for (int i = 0; i < N_bits; ++i) {
        rate10_next[i] = fr_current.rate10[i] + std_normal(rng) * step_size_schedule[step];
        rate01_next[i] = fr_current.rate01[i] + std_normal(rng) * step_size_schedule[step];
      }
      Vdbl corr1_next(corr_free), corr0_next(corr_free);
      for (int k = 0; k < corr_free; ++k) {
        double sd = step_size_schedule[step] * corr_step_scale[k];
        corr1_next[k] = fr_current.corr1[k] + std_normal(rng) * sd;
        corr0_next[k] = fr_current.corr0[k] + std_normal(rng) * sd;
      }

      // Stochastic pruning: retain current values at a randomly chosen subset of bits/pairs
      Vint retained_bits  = dg_sample_without_replacement(N_bits, n_retained_bits[step], rng);
      Vint retained_corrs = dg_sample_without_replacement(corr_free, n_retained_corrs[step], rng);
      for (int i : retained_bits) {rate10_next[i] = fr_current.rate10[i]; rate01_next[i] = fr_current.rate01[i];}
      for (int k : retained_corrs) {corr1_next[k] = fr_current.corr1[k]; corr0_next[k] = fr_current.corr0[k];}

      bool in_burnin = step > FR_search && step <= FR_search + burnin;
      bool have_means = !std::isnan(rate10_mean) && !std::isnan(rate01_mean);

      if (in_burnin && have_means) {
        double m10 = std::accumulate(rate10_next.begin(), rate10_next.end(), 0.0) / N_bits;
        double m01 = std::accumulate(rate01_next.begin(), rate01_next.end(), 0.0) / N_bits;
        for (int i = 0; i < N_bits; ++i) {
          rate10_next[i] += rate10_mean - m10;
          rate01_next[i] += rate01_mean - m01;
        }
      }
      for (int i = 0; i < N_bits; ++i) {
        bool oob10 = rate10_next[i] < 0.0 || rate10_next[i] > 1.0;
        bool oob01 = rate01_next[i] < 0.0 || rate01_next[i] > 1.0;
        if (in_burnin && have_means) {
          if (oob10) {rate10_next[i] = rate10_mean;}
          if (oob01) {rate01_next[i] = rate01_mean;}
        } else {
          if (oob10) {rate10_next[i] = fr_current.rate10[i];}
          if (oob01) {rate01_next[i] = fr_current.rate01[i];}
        }
      }
      for (int k = 0; k < corr_free; ++k) {
        if (corr1_next[k] < -1.0 || corr1_next[k] > 1.0) {corr1_next[k] = fr_current.corr1[k];}
        if (corr0_next[k] < -1.0 || corr0_next[k] > 1.0) {corr0_next[k] = fr_current.corr0[k];}
      }

      Vdbl FR_next(n_params);
      for (int i = 0; i < N_bits; ++i) {FR_next[i] = rate10_next[i]; FR_next[N_bits + i] = rate01_next[i];}
      for (int k = 0; k < corr_free; ++k) {FR_next[2*N_bits + k] = corr1_next[k]; FR_next[2*N_bits + corr_free + k] = corr0_next[k];}
      FlipRates fr_next = pack_fr(FR_next, N_bits);

      DGStepResult result = dg_simulation_step(
        STdata, fr_next, N_cells, blank_weight_schedule[step], maintain_gene_identity, (int)fresh_seed(), n_forks, rng
      );
      double sim_mse_next = result.sim_mse;
      if (w.prior > 0.0 || w.dispersion > 0.0) {
        sim_mse_next += fr_prior_penalty(fr_next, fr_priors, w.prior, w.dispersion, false, nullptr);
      }
      if (w.feasibility > 0.0) {
        sim_mse_next += w.feasibility * dg_feasibility_penalty(fr_next, STdata.cb, max_fs_bc);
      }

      double acceptance_prob = std::min(1.0, std::exp(-(sim_mse_next - sim_mse_current) / temperature_schedule[step]));
      double draw = boost::random::uniform_01<double>()(rng);
      bool   forced_accept = max_calls_per_step > 0 && attempts_this_step >= max_calls_per_step;
      if (draw < acceptance_prob || step == FR_search || forced_accept) {

        fr_current      = fr_next;
        sim_mse_current = sim_mse_next;
        step++;
        attempts_this_step = 0;

        msle[step - 1]             = sim_mse_current;
        stip_rate10_hist[step - 1] = std::accumulate(fr_next.rate10.begin(), fr_next.rate10.end(), 0.0) / N_bits;
        stip_rate01_hist[step - 1] = std::accumulate(fr_next.rate01.begin(), fr_next.rate01.end(), 0.0) / N_bits;
        for (int g = 0; g < N_genes; ++g) {PPV_genes_mat(step - 1, g) = result.PPV_genes[g];}
        if (step >= resample_window_start) {
          resample_read.push_back(result.read_counts);
          resample_corrected.push_back(result.corrected_counts);
          resample_hit.push_back(result.hit_counts);
          resample_rate10.push_back(fr_next.rate10);
          resample_rate01.push_back(fr_next.rate01);
          resample_corr1.push_back(fr_next.corr1);
          resample_corr0.push_back(fr_next.corr0);
        }

        if (step == FR_search) {
          maintain_gene_identity = true;
          int lo = std::max(1, FR_sample_start), hi = FR_search;
          double sum10 = 0.0, sum01 = 0.0; int n_samp = 0;
          for (int s = lo; s <= hi; ++s) {
            double v10 = stip_rate10_hist[s - 1], v01 = stip_rate01_hist[s - 1];
            if (!std::isnan(v10) && !std::isnan(v01)) {sum10 += v10; sum01 += v01; n_samp++;}
          }
          if (n_samp > 0) {rate10_mean = sum10 / n_samp; rate01_mean = sum01 / n_samp;}
          Rcpp::Rcout << "\n\nFlip-rate search complete, mean flip rates:"
                      << "\n\t1->0: " << std::round(rate10_mean * 1000.0) / 1000.0
                      << "\n\t0->1: " << std::round(rate01_mean * 1000.0) / 1000.0 << std::endl;
          Vdbl FR_reset = draw_random_fr(STdata, fr_priors, fresh_seed());
          fr_current = pack_fr(FR_reset, N_bits);
        }
      }

      bool at_tracker = std::find(progress_tracker.begin(), progress_tracker.end(), step) != progress_tracker.end();
      if ((at_tracker || step == 10 || step == FR_search) && step != 0 && step != last_message) {
        double acceptance_rate = static_cast<double>(step - last_step) / static_cast<double>(ctr - last_ctr);
        last_step = step; last_ctr = ctr;
        if (acceptance_rate < 0.2) {
          for (int s = step; s < n_steps; ++s) {temperature_schedule[s] *= 1.1; step_size_schedule[s] *= 0.9;}
        } else if (acceptance_rate > 0.3) {
          for (int s = step; s < n_steps; ++s) {temperature_schedule[s] *= 0.9; step_size_schedule[s] *= 1.1;}
        }
        Rcpp::Rcout << "\n\nStep " << step << " / " << n_steps << " complete"
                    << "\n\tCurrent mse: " << std::round(sim_mse_current * 1000.0) / 1000.0
                    << "\n\tAcceptance rate (this batch): " << std::round(acceptance_rate * 1000.0) / 1000.0
                    << "\n\tAcceptance rate (overall): " << std::round(static_cast<double>(step) / ctr * 1000.0) / 1000.0
                    << "\n\tCalls: " << ctr << std::endl;
        last_message = step;
      }

    }

    // Collapse the resamples-window history into misread.qc()'s own erctc_plus shape: one row per
    // barcode (genes and blanks), erc/ecc/ehc as the median simulated read/corrected/hit count over
    // the window, CR/PPV from those medians, and PPV_lower/PPV_upper as the window's own 2.5%/97.5%
    // empirical quantiles of the per-step PPV ratio -- the direct analogue of misread.qc()'s delta-
    // method CI, but read straight off the chain's actual sampling distribution instead of a local
    // Gaussian approximation. Blanks get a real erc/ecc/CR (other barcodes misread into them), but
    // ehc/PPV for blanks will be ~0 by construction: dg_draw_true_counts() always draws 0 true counts
    // for blanks, so there's never a true blank spot to be a "hit" -- matching misread.qc()'s own
    // analytic model, which zeroes blanks' ground-truth count the same way.
    NumericMatrix erctc_plus(N_barcodes, 7);
    colnames(erctc_plus) = CharacterVector({"erc", "ecc", "ehc", "CR", "PPV", "PPV_lower", "PPV_upper"});
    int n_resample_steps = resample_corrected.size();
    for (int b = 0; b < N_barcodes; ++b) {
      Vdbl erc_vals(n_resample_steps), ecc_vals(n_resample_steps), ehc_vals(n_resample_steps), ppv_ratio(n_resample_steps);
      for (int s = 0; s < n_resample_steps; ++s) {
        erc_vals[s] = resample_read[s][b];
        ecc_vals[s] = resample_corrected[s][b];
        ehc_vals[s] = resample_hit[s][b];
        ppv_ratio[s] = resample_corrected[s][b] > 0.0
          ? resample_hit[s][b] / resample_corrected[s][b]
          : std::numeric_limits<double>::quiet_NaN();
      }
      double erc_med = dg_median_skipna(erc_vals);
      double ecc_med = dg_median_skipna(ecc_vals);
      double ehc_med = dg_median_skipna(ehc_vals);
      erctc_plus(b, 0) = erc_med;
      erctc_plus(b, 1) = ecc_med;
      erctc_plus(b, 2) = ehc_med;
      erctc_plus(b, 3) = ecc_med > 0.0 ? erc_med / ecc_med : NA_REAL;
      erctc_plus(b, 4) = ecc_med > 0.0 ? ehc_med / ecc_med : NA_REAL;
      erctc_plus(b, 5) = dg_quantile_skipna(ppv_ratio, 0.025);
      erctc_plus(b, 6) = dg_quantile_skipna(ppv_ratio, 0.975);
    }

    // Collapse the same resamples window into a single "fliprates" vector -- misread.qc()'s own
    // point-estimate layout (rate10 x N_bits, rate01 x N_bits, corr1 then corr0 x corr_free), but each
    // entry the posterior median over the window rather than an L-BFGS optimum. R applies
    // misread.qc()'s own naming convention to this vector, same as it does for misread.qc()'s.
    NumericVector fliprates(n_params);
    for (int i = 0; i < N_bits; ++i) {
      Vdbl r10(n_resample_steps), r01(n_resample_steps);
      for (int s = 0; s < n_resample_steps; ++s) {r10[s] = resample_rate10[s][i]; r01[s] = resample_rate01[s][i];}
      fliprates[i]          = dg_median_skipna(r10);
      fliprates[N_bits + i] = dg_median_skipna(r01);
    }
    for (int k = 0; k < corr_free; ++k) {
      Vdbl c1(n_resample_steps), c0(n_resample_steps);
      for (int s = 0; s < n_resample_steps; ++s) {c1[s] = resample_corr1[s][k]; c0[s] = resample_corr0[s][k];}
      fliprates[2*N_bits + k]            = dg_median_skipna(c1);
      fliprates[2*N_bits + corr_free + k] = dg_median_skipna(c0);
    }

    return List::create(
      _["msle"]        = msle,
      _["PPV_by_step"] = PPV_genes_mat,
      _["erctc_plus"]  = erctc_plus,
      _["fliprates"]   = fliprates,
      _["gene_idx"]    = IntegerVector(STdata.cb.genes.begin(), STdata.cb.genes.end()),
      _["blank_idx"]   = IntegerVector(STdata.cb.blanks.begin(), STdata.cb.blanks.end())
    );

    /*
     * Literal correlated-Gaussian-noise spot simulation (see the DG correlation-matrix helpers
     * above, and make_SpotSim_DG()) -- schedules, two-phase FR_search-then-burnin design,
     * stochastic bit/pair pruning, and burn-in mean-rate enforcement, operating on FlipRates
     * (rate10/rate01/corr1/corr0) and FlipRatePriors. See misread.qc.DG()'s own documentation
     * (R/main.R) for the full argument/return contract; this function expects burnin/FR_search/
     * weights$max_blank/weights$min_blank already resolved (no NULL defaults at this level -- see
     * pack_weights()), and is not itself documented via roxygen (no man page), matching mQC()'s
     * own convention of leaving the internal C++-exported engine undocumented at the Rd level.
     *
     * max_calls_per_step (see misread.qc.DG()'s own @param doc for the user-facing explanation):
     * 0 means unlimited, matching all prior behavior exactly. When > 0, a step that fails to
     * accept a proposal within this many attempts force-accepts its last proposal and moves on
     * regardless of acceptance_prob -- bounding worst-case runtime per step, since the annealing's
     * acceptance rate can otherwise collapse for some (seed, data, schedule) combinations and
     * need many proposals to advance (each one potentially expensive -- see make_SpotSim_DG()'s
     * own neighbors). An opt-in safety valve for development/test use (this package's own fast
     * smoke tests use it), not routine production fits.
     *
     * seed: seeds the single boost::random::mt19937 (`rng`, declared just above) that drives this
     * whole walk -- the Gaussian proposal draws, the stochastic bit/pair pruning
     * (dg_sample_without_replacement()), the acceptance draw, and (via dg_draw_true_counts())
     * per-cell Gamma/Poisson transcript-count draws, plus every per-step `fresh_seed()` handed to
     * make_SpotSim_DG()'s forked children. R's own RNG/set.seed() is not touched anywhere in this
     * function, matching draw_random_fr()'s and misread.qc()'s own explicit-seed convention; a
     * caller's set.seed() in R has no effect on mQC_DG()'s/misread.qc.DG()'s reproducibility.
     */
  }

/*
 * *********************************************************************************************************************
 * Test helpers (exported for use in testthat only)
 */

// Sum TR(bc, flips, fr) over all 2^N_bits flip patterns; must equal 1.0 for every valid parameter set.
// [[Rcpp::export]]
double tr_sum_check(
    int                        bc,
    const std::vector<double>& rate10,
    const std::vector<double>& rate01,
    const std::vector<double>& corr1,
    const std::vector<double>& corr0
  ) {
    int N_bits = rate10.size();
    Vdbl params;
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
   
    /*
     * Arguments mirror pack_fr: rate10 and rate01 are per-bit flip rates; corr1 and corr0 are the
     * strict lower-triangle correlation vectors (length N_bits*(N_bits-1)/2 each). Signature
     * spelled out as std::vector<double> rather than the Vdbl typedef because
     * Rcpp::compileAttributes()/RcppExports.cpp can't resolve local typedefs from main.cpp.
     *
     * Initial draft of this function written by Claude Sonnet 4.6.
     */
  }

// Test export wrapping bvn_cdf_uppertail() for verifying against R's mvtnorm::pmvnorm
// [[Rcpp::export]]
double bvn_cdf_uppertail_test(double h, double k, double rho) {return bvn_cdf_uppertail(h, k, rho);}

// Test export wrapping dg_find_sigma_RootBisection(); see bvn_cdf_uppertail_test()'s own note.
// [[Rcpp::export]]
double dg_find_sigma_RootBisection_test(double threshold_i, double threshold_j, double desired_cov) {return dg_find_sigma_RootBisection(threshold_i, threshold_j, desired_cov);}

// Test export wrapping dg_nearPD_native(), checked against R's Matrix::nearPD() over random matrices.
// [[Rcpp::export]]
Eigen::MatrixXd dg_nearPD_native_test(Eigen::MatrixXd x) {return dg_nearPD_native(x, true, true, true);}

// Test export pinning boost::random::mt19937's exact draw sequence for a fixed seed.
// [[Rcpp::export]]
std::vector<unsigned int> rng_engine_draw_test(
    unsigned int seed, 
    int n
  ) {
    boost::random::mt19937     rng(seed);
    std::vector<unsigned int>  out(n);
    for (int i = 0; i < n; ++i) {out[i] = rng();}
    return out;
   
    /*
     * This and the other rng_*_draw_test() exports below pin the exact draw sequence
     * boost::random::mt19937 (and its normal_distribution/bernoulli_distribution) produce for a
     * fixed seed. The statistical properties of these draws are not in question; the point is to
     * catch anything that silently changes the *specific numbers* a fixed seed produces -- a
     * reversion to std::, a Boost upgrade that changes a distribution's internal sampling
     * algorithm, or a change to call order/seeding -- since any of those would make simulation
     * results produced before vs. after the change diverge for users who rely on a fixed seed for
     * reproducibility (see draw_random_fr() above).
     */
  }

// Test export pinning boost::random::normal_distribution's exact draw sequence for a fixed seed.
// [[Rcpp::export]]
std::vector<double> rng_normal_draw_test(
    unsigned int seed, 
    int n, 
    double mean, 
    double sd
  ) {
    boost::random::mt19937                     rng(seed);
    boost::random::normal_distribution<double> dist(mean, sd);
    std::vector<double>                        out(n);
    for (int i = 0; i < n; ++i) {out[i] = dist(rng);}
    return out;
  }

// Test export pinning boost::random::bernoulli_distribution's exact draw sequence for a fixed seed.
// [[Rcpp::export]]
std::vector<int> rng_bernoulli_draw_test(
    unsigned int seed, 
    int n, 
    double p
  ) {
    boost::random::mt19937                        rng(seed);
    boost::random::bernoulli_distribution<double> dist(p);
    std::vector<int>                              out(n);
    for (int i = 0; i < n; ++i) {out[i] = dist(rng) ? 1 : 0;}
    return out;
  }

// Test export pinning boost::random::gamma_distribution's exact draw sequence for a fixed seed.
// [[Rcpp::export]]
std::vector<double> rng_gamma_draw_test(
    unsigned int seed,
    int          n,
    double       shape,
    double       scale
  ) {
    boost::random::mt19937                    rng(seed);
    boost::random::gamma_distribution<double> dist(shape, scale);
    std::vector<double>                       out(n);
    for (int i = 0; i < n; ++i) {out[i] = dist(rng);}
    return out;
  }

// Test export pinning boost::random::poisson_distribution's exact draw sequence for a fixed seed.
// [[Rcpp::export]]
std::vector<int> rng_poisson_draw_test(
    unsigned int seed,
    int          n,
    double       mean
  ) {
    boost::random::mt19937                 rng(seed);
    boost::random::poisson_distribution<int> dist(mean);
    std::vector<int>                       out(n);
    for (int i = 0; i < n; ++i) {out[i] = dist(rng);}
    return out;
  }

// Test export pinning the param_type-reset reuse pattern in draw_random_fr()'s normal_distribution.
// [[Rcpp::export]]
std::vector<double> rng_normal_paramtype_draw_test(
    unsigned int                seed,
    double                      mean1, double sd1,
    double                      mean2, double sd2
  ) {
    boost::random::mt19937                     rng(seed);
    boost::random::normal_distribution<double> dist(0.0, 1.0);
    boost::random::normal_distribution<double>::param_type p1(mean1, sd1);
    boost::random::normal_distribution<double>::param_type p2(mean2, sd2);
    std::vector<double> out(2);
    out[0] = dist(rng, p1);
    out[1] = dist(rng, p2);
    return out;
   
    /*
     * One normal_distribution object, re-parameterized per draw via param_type rather than
     * reconstructed -- a distinct boost::random code path from a freshly constructed distribution,
     * worth pinning separately.
     */
  }

// Test export wrapping dg_ppois_lower() (this file's native replacement for R::ppois(x, lambda, 1,
// 0)), for direct comparison against R's own ppois() over a range of x/lambda.
// [[Rcpp::export]]
double dg_ppois_lower_test(double x, double lambda) {return dg_ppois_lower(x, lambda);}

