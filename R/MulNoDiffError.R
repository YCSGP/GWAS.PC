#' This function calculates the non-centrality parameter of the chisquare distribution
#' when the genotyping error is the same in cases and controls, 
#' assuming multiplicative mode of inheritance
#' @param r relative risk
#' @param pA minor allele frequency
#' @param ncase number of cases 
#' @param nctrl number of controls
#' @param pD disease prevalence
#' @param error.mat 3*3 matrix of genotyping error
#' @return the non-centrality parameter of the chisquare distribution 
#' @return assuming the genotyping error is the same in cases and controls
#' @export
mul.ndiff<-function(r,pA, ncase, nctrl, pD, error.mat)
{
  pa  = 1-pA
  paa = pa^2 # P(aa)
  pAa = 2*pa*pA
  pAA = pA^2
  pD.aa = pD / (paa + r*pAa + r^2*pAA) #P(D|aa)
  pD.Aa = r*pD.aa
  pD.AA = r^2*pD.aa
  pX.D = c(pD.AA, pD.Aa, pD.aa)*c(pAA, pAa, paa)/pD
  pX.Dbar = (1-c(pD.AA, pD.Aa, pD.aa))*c(pAA, pAa, paa) / (1-pD)
  pX.D=pX.D%*%error.mat
  pX.Dbar=pX.Dbar%*%error.mat
  x=c(2, 1, 0)
  y=matrix(c(pX.D*ncase, pX.Dbar*nctrl), nrow = 3, ncol = 2, byrow = FALSE)
  y<-round(y)
  fit = glm(y~x, family=binomial(link = "logit")) #Null Deviance is the noncentral param.
  fit$null.deviance
}

