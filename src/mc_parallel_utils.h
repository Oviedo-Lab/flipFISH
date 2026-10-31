#ifndef FLIPFISH_MC_PARALLEL_UTILS_H
#define FLIPFISH_MC_PARALLEL_UTILS_H

#include <Rcpp.h>
#include <functional>

/*
 * Batch map-reduce dispatched through R's own parallel::mcparallel()/mccollect() (POSIX fork(),
 * via the parallel package's own already-CRAN-accepted compiled code) rather than flipFISH
 * calling fork()/_exit() directly. Used by expected_bc_counts(), make_SpotSim(), and
 * make_SpotSim_DG() to distribute per-barcode work across processes when n_forks > 1.
 *
 * build_call(i) returns the Rcpp::Language call for batch i -- naming one of this package's own
 * "compute one batch" [[Rcpp::export]] functions plus that batch's own arguments (a shared,
 * already-built, read-only context passed as a non-owning XPtr, plus whatever per-call data --
 * the current flip rates, a random seed, etc. -- changes from one call to the next). Each caller
 * builds its own call since each use site's argument shape differs. reduce(i, result) folds
 * batch i's returned R object (whatever build_call(i)'s target function returns -- typically an
 * Rcpp::List) into the caller's own output.
 *
 * Callers must only invoke this with n_use > 1 on POSIX platforms: parallel::mcparallel() has no
 * Windows equivalent (no fork()), so every call site guards n_use to 1 there (the same `#ifdef
 * _WIN32` pattern used throughout this package) before it ever reaches here.
 */
inline void run_mc_batches(
    int                                 n_use,
    std::function<Rcpp::Language(int)>  build_call,
    std::function<void(int, SEXP)>      reduce
  ) {
    Rcpp::Environment parallel_ns = Rcpp::Environment::namespace_env("parallel");
    Rcpp::Function    mcparallel  = parallel_ns["mcparallel"];
    Rcpp::Function    mccollect   = parallel_ns["mccollect"];

    Rcpp::List jobs(n_use);
    for (int i = 0; i < n_use; ++i) {
      jobs[i] = mcparallel(build_call(i));
    }
    Rcpp::List results = mccollect(jobs);
    for (int i = 0; i < n_use; ++i) {
      reduce(i, results[i]);
    }
  }

#endif // FLIPFISH_MC_PARALLEL_UTILS_H
