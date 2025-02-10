
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCAssess

<!-- badges: start -->
<!-- badges: end -->

PCAssess is an R package for assessing high grading bias in population
genetics. PCAssess implements and automates permuation tests to detect
high grading bias in Principal Component Analysis (PCAs).

## Installation

You can install the development version of PCAssess from
[GitHub](https://github.com/hemstrow/PCAssess) with:

``` r
# install.packages("pak")
pak::pak("hemstrow/PCAssess")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PCAssess)
#> Loading required package: data.table
#> Warning: package 'data.table' was built under R version 4.3.3
#> Loading required package: foreach
## basic example code
#perm_res <- run_permutation(mon_sn, "pop", 10, fst_cut=.95, par = FALSE, store_pca = FALSE) # run permutation test on example Monarch data
#plot_permutation_res(perm_res, 10) # plot observed PCA and 100 permutations
#plot_permutation_res(perm_res, 0) # plot only the observed PCA
```
