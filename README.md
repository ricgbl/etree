
<!-- README.md is generated from README.Rmd. Please edit that file -->

# etree

<!-- badges: start -->
<!-- badges: end -->

The goal of etree is to provide a friendly implementation of Energy
Trees, a model for classification and regression with structured and
mixed-type data. The package currently cover functions and graphs as
structured covariates.

## Installation

You can install the development version of etree from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ricgbl/etree")
```

## Example

This is a basic example which shows how to fit an Energy Tree for
regression using a toy dataset with four covariates of different types:
numeric, nominal, functional, and in the form of graphs.

``` r
library(etree)

# Covariates
nobs <- 100
cov_num <- rnorm(nobs)
cov_nom <- factor(rbinom(nobs, size = 1, prob = 0.5))
cov_gph <- lapply(1:nobs, function(j) igraph::sample_gnp(100, 0.2))
cov_fun <- fda.usc::rproc2fdata(nobs, seq(0, 1, len = 100), sigma = 1)
cov_list <- list(cov_num, cov_nom, cov_gph, cov_fun)
 
# Response variable
resp_reg <- cov_num ^ 2

# Energy Tree fit
etree_fit <- etree(response = resp_reg, 
                   covariates = cov_list)
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> Warning in .create_newcov(covariates = covariates, response = response, : No
#> names available for covariates. Numbers are used instead.
```

Additional and more complex examples can be found in the package’s
vignettes.
