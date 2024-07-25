# scNMF

<!-- badges: start -->

<!-- badges: end -->

The goal of scNMF is to run NMF programs analysis on single cell sequencing data ...

## Installation

You can install the development version of scNMF like so:

``` r
remotes::install_github('Moonerss/scNMF')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scNMF)

## a list of cpm matrix: `matrix_list`
nmf_programs_list <- nmf_programs(matrix_list, rank = 2:9)

top_genes_list <- select_top_program_genes(nmf_programs_list, top_n = 50)

robust_nmf_list <- select_robust_nmf_program(top_genes_list)

p <- plot_similarity_matrix(top_genes_list, robust_nmf_list)

```
