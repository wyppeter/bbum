
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bbum - BBUM correction for significance testing of *p* values

## Description

Standard multiple testing correction methods cannot distinguish dominant
primary effects from background secondary effects above null. A
bi-beta-uniform mixture (BBUM) model allows the modeling and correction
for the false discovery rate (FDR) to exclude both null and secondary
effects.

------------------------------------------------------------------------

## Installation

You can install the `bbum` package from
[GitHub](https://github.com/wyppeter/bbum) with:

``` r
if(!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("wyppeter/bbum")
library(bbum)
```

------------------------------------------------------------------------

## Functions

The `bbum` package defines standard `d`, `p`, and `r` functions for
bi-beta-uniform mixture (BBUM) models as well as the original
beta-uniform mixture (BUM) models. `q` functions (quantile functions)
cannot be easily expressed for BUMs and BBUMs and are unavailable.

``` r
# BBUM distribution ----
dbbum(x = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02, theta = 0.11, r = 0.003)
pbbum(q = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02, theta = 0.11, r = 0.003)
## Random generator ----
rbbum(n = 100, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
# This generates values with their associated distribution components
rbbum.ID(n = 100, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)

# BUM distribution ----
dbum(x = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02)
pbum(q = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02)
```

Typical usage for p value fitting and/or correction would make use of
one of these three functions, depending on the level of work needed:

``` r
# BBUM_fit: Simple fitting ----
# Used similarly to, say, lm()
BBUM_fit(
  dt_signal_set = c(0.000021, 0.00010, 0.03910, 0.031, 0.001,
                    0.13, 0.21, 0.38, 0.42, 0.52, 0.60, 0.73, 0.81, 0.97),
  dt_bg_set     = c(0.501, 0.203, 0.109, 0.071, 0.019,
                    0.11, 0.27, 0.36, 0.43, 0.50, 0.61, 0.77, 0.87, 0.91),
  starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1))
)
```

``` r
# BBUM_corr: p value adjustment/correction for FDR ----
# Used similarly to, say, p.adjust()
BBUM_corr(
  pvals         = c(0.501, 0.203, 0.109, 0.071, 0.019, 0.031, 0.001,
                    0.000021, 0.00010, 0.03910,
                    0.11, 0.27, 0.36, 0.43, 0.50, 0.61, 0.77, 0.87, 0.91,
                    0.13, 0.21, 0.38, 0.42, 0.52, 0.60, 0.73, 0.81, 0.97),
  signal_set    = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
                    TRUE,  TRUE,  TRUE,
                    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                    TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE)
)
```

``` r
# BBUM_DEcorr: DESeq2 table processing and significance calling ----
# Used similarly to, say, lfcShrink()
# DESeq2 calls
dds = DESeqDataSetFromMatrix(countData = cts, ...)
dds = DESeq(dds)
res = results(dds) %>%
  as.data.frame() %>%
  mutate(FCup = log2FoldChange > 0)
# Then call BBUM_DEcorr
res.BBUMcorr = BBUM_DEcorr(
  df.deseq = res,
  classCol = "FCup"
  )
```

The package also includes many `ggplot2` plots for QC and analysis of
`BBUM_DEcorr` results:

``` r
# BBUM_plot: DE results plotting ----
BBUM_plot(df.bbum = res.BBUMcorr,
          option = "hist"
          )
```

For more details on included functions and how to use them, refer to
documentation in the package:

``` r
# For example
?pbbum
?BBUM_fit
?BBUM_corr
?BBUM_DEcorr
```

------------------------------------------------------------------------

## References

This package was created by Peter Y. Wang (@wyppeter) in David Bartel’s
lab. For more details on the theoretical basis, usage, and benchmarking
of the package, refer to the preprint at:

> Wang PY, Bartel DP. 2021. TBD. *bioRxiv*. <doi:TBD>
