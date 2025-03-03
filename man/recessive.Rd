% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RecessiveNoError.R
\name{recessive}
\alias{recessive}
\title{This function returns the non-centrality parameter when gentoyping error is ignored,
assuming recessive mode of inheritance}
\usage{
recessive(r, pA, ncase, nctrl, pD)
}
\arguments{
\item{r}{relative risk}

\item{pA}{minor allele frequency}

\item{ncase}{number of cases}

\item{nctrl}{number of control samples}

\item{pD}{disease prevalence}
}
\value{
non-centrality parameter of the chi-square distribution under alternative hypothesis
}
\description{
This function returns the non-centrality parameter when gentoyping error is ignored,
assuming recessive mode of inheritance
}

