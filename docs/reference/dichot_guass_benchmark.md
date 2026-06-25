# Benchmark misread QC with dichotomized Gaussian simulation

This function takes the same summary statistics (`STdata`) and barcode
codebook (`codebook`) as `misread.qc` and runs dichotomized Gaussian
simulations with stipulated bit-flip rates and bit-flip correlations in
order to estimate how well the L-BFGS algorithm recovers the bit-flip
rates and bit-flip correlations for the given data set and codebook.

## Usage

``` r
dichot_guass_benchmark(
  STdata,
  codebook,
  n_sims = 100,
  n_forks = 1,
  max_flips = 0,
  report_freq = 10,
  maxeval = 500,
  max_correctable_Hamming_distance = NULL
)
```

## Arguments

- STdata:

  Numeric matrix with rows as barcodes, columns labeled "rates",
  "variance", "counts". Must have barcode names (e.g., gene or protein
  species) as row names.

- codebook:

  Codebook with barcodes as row names and bits as columns.

- n_sims:

  Number of dichotomized Gaussian simulations to run. The default is
  100.

- n_forks:

  Number of parallel forks to use for expected-count computation,
  default is 1. Must be 1 on Windows, can be higher on Mac and Linux.

- max_flips:

  When analytically computing expected corrected counts per barcode, the
  function will ignore misreads larger than this hamming distance,
  default is 0, which is interpreted as no limit. Using all misreads
  will likely be prohibitively slow; a value between six and ten is
  probably advisable. Values of 3 or 4 work well for initial trouble
  shooting and testing.

- report_freq:

  Divisor specifying report frequency during optimization; will print
  updates every report_freq accepted calls, default 10.

- maxeval:

  Maximum number of objective function evaluations for L-BFGS, default
  500.

- max_correctable_Hamming_distance:

  Maximum Hamming distance for misreads to be corrected, default NULL
  sets it to one less than the minimum Hamming distance between codebook
  entries.

## Value

A list giving::

- `fr_est`: A matrix giving the estimated flip rates and bit-flip
  correlations (columns) for each simulation run (rows).

- `fr_stipulated`: A vector giving the stipulated flip rates and
  bit-flip correlations used to generate the dichotomized Gaussian
  simulations.

- `PPV_est`: A matrix giving positive predictive value (PPV) per barcode
  (columns), estimated by L-BFGS each simulation (rows).

- `PPV_expected`: A vector giving the PPV values expected via analytical
  computation, given the stipulated flip rates and stipulated bit-flip
  correlations in `fr_stipulated`.

- `ecc_est`: A matrix giving the expected corrected count per barcode
  (columns), estimated by L-BFGS each simulation (rows).

- `ecc_expected`: A vector giving the corrected count per barcode
  expected via analytical computation, given the stipulated flip rates
  and stipulated bit-flip correlations in `fr_stipulated`.

- `sim_counts`: A matrix giving the actual simulated count per barcode
  (columns) for each simulation (rows).
