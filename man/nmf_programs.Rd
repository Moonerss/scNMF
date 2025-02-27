% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{nmf_programs}
\alias{nmf_programs}
\title{Run NMF}
\usage{
nmf_programs(
  mat,
  lognormalize = T,
  do_filtering = TRUE,
  min_expr = 3.5,
  min_pct = 0.2,
  center = TRUE,
  scale = FALSE,
  non_negative = TRUE,
  rank = 2,
  method = "snmf/r",
  seed = 1,
  ...
)
}
\arguments{
\item{mat}{a numeric matrix with genes in row and cells in column.}

\item{lognormalize}{whether log transformed the data, if TRUE, the data will transformed with \code{log2(x+1)}}

\item{do_filtering}{Whether filter the gene by \code{min_expr} and \code{min_pct}.}

\item{min_expr}{The minum expression of gene expression}

\item{min_pct}{The minum precent of cell expression gene greater than \code{min_expr}.}

\item{center}{Whether to center the data matrix}

\item{scale}{Whether to scale the data matrix}

\item{non_negative}{Whether fill the negative value with zero}

\item{rank}{specification of the factorization rank. Here it must be a single numeric value.}

\item{method}{specification of the NMF algorithm.}

\item{seed}{specification of the starting point or seeding method.}

\item{...}{Other argument in \code{nmf} function.}
}
\value{
return a list contain w and h matrix of NMF.
}
\description{
Run nonnegative matrix factorization (NMF) with matrix
}
