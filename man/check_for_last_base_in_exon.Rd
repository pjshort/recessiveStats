% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/tweak_lofs.R
\name{check_for_last_base_in_exon}
\alias{check_for_last_base_in_exon}
\title{modifies the consequence where the variant is at the last base in the exon}
\usage{
check_for_last_base_in_exon(variant, exon_ends, strand)
}
\arguments{
\item{variant}{row of dataframe}

\item{exon_ends}{vector of exon end positions for the variants gene.}

\item{strand}{which chromosome strand the gene is on, 1 for plus strand, and
-1 for minus strand.}
}
\value{
consequence for variant.
}
\description{
modifies the consequence where the variant is at the last base in the exon
}

