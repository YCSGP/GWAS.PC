% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiNoDiffError.R
\name{NCP}
\alias{NCP}
\title{This function calculates the non-centrality parameter of burden test statistic under alternative hypothesis
when the genotyping error is the same in cases and controls.}
\usage{
NCP(r, pA, pD, m10, m00, M, Ncase)
}
\arguments{
\item{r}{relative risk}

\item{pA}{minor allele frequency}

\item{pD}{disease prevalence}

\item{m10}{1-sensitivity}

\item{m00}{specificity}

\item{M}{number of markers in a gene/region}

\item{Ncase}{number of cases}
}
\value{
the non-centrality parameter of burden test statistic under alternative hypothesis
}
\description{
This function calculates the non-centrality parameter of burden test statistic under alternative hypothesis
when the genotyping error is the same in cases and controls.
}

