# multibiasmeta

<!-- badges: start -->
[![R-CMD-check](https://github.com/mikabr/multibiasmeta/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mikabr/multibiasmeta/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`multibiasmeta` is an R package that provides bias correction and sensitivity analysis for the joint effects of within-study and across-study biases in meta-analysis (per [Mathur, 2022](doi:10.31219/osf.io/u7vcb)).


## Installation

You can install the development version of multibiasmeta from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mikabr/multibiasmeta")
```

## Example

Calculate the meta-analytic effect size estimate, correcting for 1) publication
bias with an assumed selection ratio of 4 (i.e., affirmative results are 4x more
likely to be published than nonaffirmative ones); 2) internal biases with
assumed mean values (on the same scale as `yi` values).

``` r
library(multibiasmeta)
multibias_meta(yi = meta_meat$yi, vi = meta_meat$vi, biased = TRUE, selection_ratio = 4,
               bias_affirmative = log(1.5), bias_nonaffirmative = log(1.1))
```

Calculate how high mean internal bias would need to be to attenuate the effect
size estimate to the null, assuming a selection ratio of 4.

``` r
multibias_evalue(yi = meta_meat$yi, vi = meta_meat$vi, selection_ratio = 4,
                 biased = !meta_meat$randomized)
```
