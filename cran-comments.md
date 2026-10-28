# CRAN comments

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.2
* (win-builder / R-hub / Debian/Linux CI not yet run -- to be added before actual submission)

## R CMD check results

`R CMD check --as-cran` currently reports 3 WARNINGs and 6 NOTEs. Notes on each:

### WARNING: compiled code calls `_exit`

```
File 'flipFISH/libs/flipFISH.so':
  Found '__exit', possibly from '_exit' (C), '_Exit' (C)
    Object: 'main.o'
```

This comes from `src/fork_utils.h`'s `run_forked_batches()`, a small fork()/pipe() map-reduce
helper used by the handful of `src/main.cpp` functions that distribute per-barcode work across
processes (`expected_bc_counts()`, `make_SpotSim()`, `make_SpotSim_DG()`). Each of these is called
on every evaluation of the package's optimizers (every L-BFGS objective evaluation in
`misread.qc()`, every MCMCSA step in `misread.qc.DG()`), so this is a hot path, not an incidental
one.

**Why `_exit()` specifically, and why it is safe here.** After `fork()`, each child computes its
batch and writes the result to its pipe; `_exit()` (not `exit()`) is used deliberately so the
child skips C++ static destructors, `atexit()` handlers, and C stdio buffer flushes it shares with
the parent -- running any of those twice (once implicitly via the parent's own eventual exit) is
the actual hazard `_exit()` avoids. The child never calls back into R: every `Rcpp::stop()`/R API
call in the file is confined to the parent process (verified by inspection -- the only two calls
inside the `if (pid == 0)` child branch are `write()` and `_exit()` themselves). Windows is handled
by `#ifdef _WIN32` excluding the whole fork()/pipe() path from compilation there -- every batch
just runs serially in-process on Windows, which `R/main.R`'s `check_forks()` also enforces at the R
level (`n_forks` is forced to 1 off of POSIX).

**Why not use `parallel::mclapply()` instead.** This was evaluated directly rather than assumed.
Benchmarking `run_forked_batches()` against an equivalent `parallel::mclapply()` implementation,
batching the same payload sizes:

* On realistic per-call workloads (the package's own 819-barcode/28-bit reference codebook), the
  hand-rolled path is 5-30% faster than the `mclapply()` equivalent at matched fork counts, because
  it moves raw binary buffers through pipes instead of paying R's `serialize()`/`unserialize()`
  cost per worker.
* On small per-call workloads -- common here, since both optimizers call this on every single
  iteration -- `mclapply()`'s per-worker overhead dominates and can make parallelizing actively
  counterproductive (observed: `mclapply()` at 8 workers *slower* than running serially, 0.0071s vs
  0.0024s, while the hand-rolled path stayed flat at ~0.0019s regardless of fork count).

Given `expected_bc_counts()` alone runs on every optimizer iteration (hundreds of times per fit),
these per-call differences compound across a full fit. The custom mechanism is a deliberate,
measured choice for this access pattern, not premature optimization, and the compiled-code WARNING
above is understood and accepted as the tradeoff for it.

### WARNING: no examples, no tests, no vignettes run during check

A `testthat` suite exists (`tests/testthat/`, 5 files) but is not yet wired up to run under
`R CMD check` (missing `tests/testthat.R`); vignettes were deliberately moved to
`vignettes/articles/` (pkgdown-only, not built as package vignettes, since their content is still
under active revision) and no `@examples` have been added yet. Known outstanding work, not yet
resolved.

### NOTE: no `Authors@R`, Title not in title case, `License` field non-standard

Known outstanding work as of this check run; being addressed.

### NOTE: `nloptr` listed in Imports but not imported from

```
Namespaces in Imports field not imported from:
  'MASS' 'Matrix' 'dplyr' 'grid' 'gridExtra' 'nloptr'
```

Expected and intentional for `nloptr`: the package's C++ code calls NLopt only through nloptr's
registered C routines (`nloptrAPI.h`, resolved via `R_GetCCallable()` at first use) rather than
linking a system `libnlopt`, so there is no R-level symbol to `importFrom`. `nloptr` is listed in
`Imports` (and force-loaded in `.onLoad()`) purely so its namespace -- and therefore its C-callable
registrations -- is guaranteed loaded before any compiled code needs it. `MASS`, `Matrix`,
`dplyr`, `grid`, `gridExtra` are leftover/unused and will be removed or wired up in a follow-up
pass.

### NOTE: `LICENSE` file / license specification

Known outstanding work as of this check run; being addressed.

### NOTE/INFO: `RcppEigen` in Imports should probably be LinkingTo-only; HTML tidy version

Pre-existing, not yet addressed (`RcppEigen`); the HTML Tidy NOTE is a limitation of the local
check environment, not the package.

## Downstream dependencies

None (new submission).
