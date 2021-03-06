% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/alleles_recessive_test.R
\name{test_enrichment}
\alias{test_enrichment}
\title{test for enrichment of inherited variants}
\usage{
test_enrichment(freq, biallelic_lof, biallelic_func, lof_func, cohort_n,
  autozygous_rate = 0)
}
\arguments{
\item{freq}{list of cumulative frequencies of variation in a population for
rare LoF variants, and rare functional variants.}

\item{biallelic_lof}{number of probands with inherited Lof/LoF variants in
the gene.}

\item{biallelic_func}{number of probands with inherited func/func variants in
the gene.}

\item{lof_func}{number of probands with inherited Lof/Func variants in the gene.}

\item{cohort_n}{number of probands in population.}

\item{autozygous_rate}{rate of autozygosity in the cohort being investigated.}
}
\value{
a list of P values from tests, under LoF and functional tests.
}
\description{
test for enrichment of inherited variants
}

