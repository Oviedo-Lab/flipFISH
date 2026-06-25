# Run misread QC analysis on summary stats data from spatial transcriptomics experiment

This function takes summary statistics from a FISH-based spatial
transcriptomics experiment and the barcode codebook (including blanks
labelled with "Blank") and runs L-BFGS (via nlopt) to make a best-fit
estimate of the bit-flip rates and correlations. The estimation (i.e.,
L-BFGS optimization) uses an analytic conditional probability model to
compute, for each barcode *b*, values expected based on bit-flip rates
and bit-flip correlations, for:

- **Read Count**: The number of spots (the "count") labelled *b* before
  error correction.

- **Corrected Count**: The number of spots (the "count") labelled *b*
  after error correction.

- **Hit Count**: The number of spots (the "count") labelled *b* after
  error correction which are in fact *b*, i.e., the number of reads
  error-corrected to *b* which are correct, aka a "hit".

The "best fit" is defined as least mean squared log error, mean squared
log error being computed by comparing the analytically implied expected
corrected counts to the observed corrected counts included in the
summary statistics. Additionally, both the expected "confidence ratio"
(CR) and expected positive predictive value (PPV) are computed, from
these implied values, for each barcode. The CR is a quality control
metric proposed by the original developers of MERFISH (e.g., see DOI
10.1016/bs.mie.2016.03.020) defined per barcode as the read count over
the corrected count, while PPV is a oft-used metric defined as the hit
count over the corrected count. That is, PPV tells us, for each barcode
*b*, the expected percentage of spots decoded as *b* which are in fact
*b*.

## Usage

``` r
misread.qc(
  STdata,
  codebook,
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

  Numeric matrix with barcodes as row names and bits as columns. All
  entries should be 1 or 0, depending on whether an mRNA molecule of the
  species represented by the row is expected to luminescence in the bit
  represented by the column. All row names from `STdata` must be
  included in the row names for `codebook`.

- n_forks:

  Number of process forks to use for expected-count computation (i.e.,
  parallel computation), default is 1. Must be 1 on Windows, can be
  higher on Mac and Linux.

- max_flips:

  When analytically computing expected corrected counts per barcode, the
  function will ignore misreads larger than this Hamming distance. The
  default is 0, which is interpreted as no limit. Using all misreads
  will likely be prohibitively slow; a value between six and ten is
  probably advisable. Values of 3 or 4 work well for initial trouble
  shooting and testing.

- report_freq:

  Divisor specifying report frequency during optimization; will print
  updates every `report_freq` accepted calls, default 10.

- maxeval:

  Maximum number of objective function evaluations for L-BFGS, default
  500.

- max_correctable_Hamming_distance:

  Maximum Hamming distance for misreads to be corrected, default NULL
  sets it to one less than the minimum Hamming distance between codebook
  entries.

## Value

A list giving:

- `STdata`: a dataframe giving the summary data from the ST run used in
  the estimation.

- `fliprates`: a labeled vector giving the estimated (i.e., best fit)
  flip rates and bit-flip correlations from the optimization.

- `erctc_plus`: a data frame giving, for each barcode, the read,
  corrected, and hit counts as well as `CR` and `PPV` values implied by
  the analytic conditional probability model, given the values in
  `fliprates`. Column names are: "erc", "ecc", "etc", "CR", and "PPV".

- `msle`: A numeric value giving the minimal mean squared log error
  found by the L-BFGS optimization of the flip rates and bit-flip
  correlations.
