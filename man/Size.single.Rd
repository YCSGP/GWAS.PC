% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MinSize.R
\name{Size.single}
\alias{Size.single}
\title{This function calculates the minimum sample size needed to achieve a specified level 
of power for single marker analysis.}
\usage{
Size.single(alpha, p, r, pA, pD, emat1, emat2 = emat1, type = "ndiff",
  moi = "m", range = c(1, 1e+05))
}
\arguments{
\item{alpha}{type I error}

\item{p}{pre-specified statistical power}

\item{r}{relative risk}

\item{pA}{minor allele frequency}

\item{pD}{disease prevalence}

\item{type}{"ndiff"/"diff", indicating non-differential or differential genotyping error}

\item{moi}{"m"/"r"/"d", indicating multiplicative, or recessive, or dominant mode of inheritance}

\item{error.mat1}{3*3 matrix of genotyping error in cases}

\item{error.mat2}{3*3 matrix of genotyping error in controls}
}
\value{
the minimum sample size needed to achieve a specified level

of power for single marker analysis.
}
\description{
This function calculates the minimum sample size needed to achieve a specified level 
of power for single marker analysis.
}

