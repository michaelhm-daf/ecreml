
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ecreml R-package

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The `ecreml` `R`-package consists of a series of functions to aid
agricultural scientists and data analysts implement the statistical
methodology described in [Mumford et
al. (2023)](https://www.sciencedirect.com/science/article/pii/S037842902300326X),
which incorporates environmental covariates (ECs) into a
multi-environment trial (MET) analysis to better understand the
environmental drivers contributing to the genotype $\times$ environment
$\times$ management practice (G $\times$ E $\times$ M) interaction. The
overall aim of the `ecreml` `R`-package is to make it as easy as
possible for scientists to implement the statistical methodology for
their own MET data.

## Installation

You can install the development version of ecreml from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("michaelhm-daf/ecreml", build_vignettes = TRUE)

```
