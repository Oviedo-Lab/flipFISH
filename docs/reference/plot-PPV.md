# Plot estimated PPV by barcode from `misread.qc` results

This function takes the results from the misread.qc function and makes a
plot of the estimated positive predictive value (PPV) for each barcode.
Barcodes are sorted by PPV value, in decreasing order, and a dashed red
line indicates the minimum PPV cutoff for "good" barcodes. Barcodes
above the cutoff are colored blue, while those that are not are colored
black.

## Usage

``` r
plot.PPV(
 qc,
 min_PPV = 0.8
)
```

## Arguments

- qc:

  List, results from misread.qc function.

- min_PPV:

  Numeric, PPV cutoff for "good" barcodes, defaults to 0.8.

## Value

ggplot object showing estimated PPV for each barcode with cutoff.
