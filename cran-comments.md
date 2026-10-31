# CRAN comments

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.2
* (win-builder / R-hub / Debian/Linux CI not yet run -- to be added before actual submission)

## R CMD check results

`R CMD check --as-cran` currently reports 1 WARNING and 5 NOTEs. Notes on each:

### WARNING: compiler warning during install (local environment artifact)

```
/Library/Frameworks/R.framework/Resources/include/R_ext/Boolean.h:62:36: warning: unknown warning
group '-Wfixed-enum-extension', ignored [-Wunknown-warning-option]
```

This originates from R's own installed header (`R_ext/Boolean.h`), not from this package's code or
`src/Makevars` (which sets no custom warning flags). It stems from a version mismatch on this local
machine between the Apple clang R itself was built with and the newer Apple clang on `$PATH` used
to compile packages (Xcode has been updated since this R binary was built); the header's
`#pragma clang diagnostic ignored "-Wfixed-enum-extension"` targets a warning-group name the newer
clang no longer recognizes under that spelling. Not expected to reproduce on CRAN's own build
machines, where the R binary and compiler toolchain are matched; will be confirmed via
win-builder/R-hub/mac-builder before submission.

### NOTE: `nloptr` listed in Imports but not imported from

```
Namespace in Imports field not imported from: 'nloptr'
```

Expected and intentional: the package's C++ code calls NLopt only through nloptr's registered C
routines (`nloptrAPI.h`, resolved via `R_GetCCallable()` at first use) rather than linking a system
`libnlopt`, so there is no R-level symbol to `importFrom`. `nloptr` is listed in `Imports` (and
force-loaded in `.onLoad()`) purely so its namespace -- and therefore its C-callable registrations
-- is guaranteed loaded before any compiled code needs it.

(`MASS`, `Matrix`, `dplyr`, `grid`, and `gridExtra` were previously also listed in `Imports` without
being used anywhere in the package's R or C++ code -- leftovers from code that no longer exists.
Removed; `Matrix` moved to `Suggests` since it is still used, but only in
`tests/testthat/test-nearPD.R`, as a reference implementation to check the package's native C++
port of `Matrix::nearPD()` against -- guarded there with `skip_if_not_installed("Matrix")`.
`RcppEigen` is likewise now `LinkingTo`-only rather than also `Imports`/`@import`, since nothing in
the package's R code calls an RcppEigen R-level function; it is used purely as C++ headers. The
`Date` field has been dropped rather than fixed, to stop it going stale again between now and
actual submission -- it will be resolved by having no stale value to notice.)

### NOTE: "New submission" / top-level files / HTML Tidy

"New submission" (under "CRAN incoming feasibility") is expected and informational for a first
submission. "Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed"
and the HTML Tidy NOTE are both limitations of the local check environment (no pandoc / no recent
HTML Tidy installed here), not the package.

### NOTE: R code for possible problems

Entirely `plot.MCMCsteps` (the unfinished function explicitly out of scope for this pass) -- its
undefined globals (`.DG_title_size`, `sim_mse`, etc.) and missing `ggplot2`/`stats` bindings will
be resolved when it is finished and exported like the other `plot_*` functions.

## Downstream dependencies

None (new submission).
