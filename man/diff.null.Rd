% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Power.R
\name{diff.null}
\alias{diff.null}
\title{This function calculates the non-centrality parameter of the null distribution of single marker test
in existence of differential genotyping error}
\usage{
\method{diff}{null}(pA, ncase, nctrl, error.mat.case, error.mat.ctrl)
}
\arguments{
\item{pA}{minor allele frequency of the SNP}

\item{ncase}{number of cases in the study}

\item{nctrl}{number of contrls in the study}

\item{error.mat.case}{3*3 matrix of genotyping error in cases}

\item{error.mat.ctrl}{3*3 matrix of genotyping error in controls}
}
\value{
non-centrality parameter of the null distribution
}
\description{
This function calculates the non-centrality parameter of the null distribution of single marker test
in existence of differential genotyping error
}

