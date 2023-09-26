
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bbum <img src='man/figures/logo.png' align="right" height="139" />

BBUM correction for significance testing of p values

## Description

Standard multiple testing correction methods are not sufficient to
directly handle datasets that contain a weaker background secondary
signal confounding the primary signal of interest. This package
describes and applies the bi-beta-uniform mixture (BBUM) model to adjust
p values. It allows a modified estimation of false discovery rate (FDR)
that takes into consideration the secondary signal in the background
that needs to be excluded.

BBUM requires two datasets or one dataset split into two subsets. One
set (*signal set*, or *sample set*) contains a mixture of the desired
primary signal points, secondary signal points, and null data points in
some distribution. The other set (*background set*) contains a mixture
of data points in an analogous distribution, except from the primary
signal. This allows the modeling of the background distribution using
the background set to assess the signal set under FDR.

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
bi-beta-uniform mixture (BBUM) model as well as the original
beta-uniform mixture (BUM) model. `q` functions (quantile functions)
cannot be easily expressed for BUM and BBUM and are unavailable.

``` r
# BBUM distribution ----
dbbum(x = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02, theta = 0.11, r = 0.003)
pbbum(q = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02, theta = 0.11, r = 0.003)
## Random generator ----
rbbum(n = 100, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)
# This generates values with their associated distribution components, as a 
# data frame with two columns
rbbum.ID(n = 100, lambda = 0.65, a = 0.1, theta = 0.02, r = 0.07)

# BUM distribution ----
dbum(x = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02)
pbum(q = c(0.013, 0.04, 0.93, 0.8), lambda = 0.73, a = 0.02)
```

Typical usage for p value fitting and/or correction would make use of
one of these three functions, depending on the level of work needed.

If you only want a BBUM model fit based on a dataset, call `BBUM_fit`.

``` r
# BBUM_fit: Simple fitting to the model ----
# Used similarly to, say, lm()
BBUM_fit(
  dt_signal_set = c(0.000021, 0.00010, 0.03910, 0.031, 0.001,
                    0.13, 0.21, 0.38, 0.42, 0.52, 0.60, 0.73, 0.81, 0.97),
  dt_bg_set     = c(0.501, 0.203, 0.109, 0.071, 0.019,
                    0.11, 0.27, 0.36, 0.43, 0.50, 0.61, 0.77, 0.87, 0.91),
  starts = list(c(lambda = 0.9, a = 0.6, theta = 0.1, r = 0.1))
)
```

If you want a BBUM-FDR-adjusted p values (pBBUM) for a set of p-values,
call `BBUM_corr`. It invokes `BBUM_fit` and then calculates the FDR, and
will be only function you need.

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

If you want BBUM-FDR correction of p-values in a DESeq dataset table
(name, log fold change, raw p-values, and a custom classifying column to
divide the data in half), call `BBUM_DEcorr`. It invokes `BBUM_corr` and
will be only function you need.

``` r
# BBUM_DEcorr: Data frame table processing and significance calling ----
# Used similarly to, say, lfcShrink(), if using DESeq2
# Process a DE results table, including p values and a user-defined custom 
#   column classifying the rows into the signal set and the background set
# DESeq2 (as example, but any similar results table should be work. See documentation.)
dds = DESeqDataSetFromMatrix(countData = cts, ...)
dds = DESeq(dds)
res = results(dds) %>%
  as.data.frame() %>%
  # For example, if only one fold change direction is possible for the primary signal
  mutate(FCup = log2FoldChange > 0)  
BBUM_DEcorr(
  df.deseq = res,
  classCol = "FCup"  # used to split dataset into two halves, one for modeling, one for analysis
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
?BBUM_plot
```

------------------------------------------------------------------------

## References

This package was created by Peter Y. Wang (@wyppeter) in David Bartel’s
lab in the Whitehead Institute. For more details on the theoretical
basis, usage, and benchmarking of the package, refer to the paper at:

> Wang PY, Bartel DP. 2023. A statistical approach for identifying
> primary substrates of ZSWIM8-mediated microRNA degradation in
> small-RNA sequencing data. *BMC Bioinformatics* **24**:195.
> <doi:10.1186/s12859-023-05306-z>

References:

> Pounds S, Morris SW. 2003. Estimating the occurrence of false
> positives and false negatives in microarray studies by approximating
> and partitioning the empirical distribution of p-values.
> *Bioinformatics* **19**:1236–1242. <doi:10.1093/bioinformatics/btg148>

> Markitsis A, Lai Y. 2010. A censored beta mixture model for the
> estimation of the proportion of non-differentially expressed genes.
> *Bioinformatics* **26**:640–646. <doi:10.1093/bioinformatics/btq001>
