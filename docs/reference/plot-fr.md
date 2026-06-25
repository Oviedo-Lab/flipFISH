# Plot estimated flip rates from `misread.qc` results

This function takes the results from the misread.qc function and makes a
plot showing the estimated bit-flip rates for each bit, separated by
flip type (1\>0 vs 0\>1).

## Usage

``` r
plot.fr(qc)
```

## Arguments

- qc:

  List, results from misread.qc function.

## Value

ggplot object showing distributions of estimated flip rates for each bit
and flip type.
