# Plot predicted vs observed counts from `misread.qc` results

This function takes the results from the misread.qc function and makes a
plot comparing the predicted error-corrected counts (ecc) to the
observed counts for each barcode. Both predicted and observed counts are
plotted on a log scale for better visibility. Barcodes are sorted by
decreasing observed count, with gene barcodes shown first and blank
barcodes shown second. A dashed vertical line indicates the separation
between gene and blank barcodes.

## Usage

``` r
plot.counts(qc)
```

## Arguments

- qc:

  List, results from misread.qc function.

## Value

ggplot object showing predicted vs observed counts for each barcode.
