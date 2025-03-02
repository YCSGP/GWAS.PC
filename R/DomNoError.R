#' This function returns the non-centrality parameter when gentoyping error is ignored,
#' assuming dominant mode of inheritance
#' @param r relative risk
#' @param pA minor allele frequency
#' @param ncase number of cases
#' @param nctrl number of control samples
#' @param pD disease prevalence
#' @return non-centrality parameter of the chi-square distribution under alternative hypothesis
#' @export
dominant<-function(r,pA, ncase, nctrl,pD)
{
  ### A is the risk allele
  pa  = 1-pA
  paa = pa^2 # P(aa)
  pAa = 2*pa*pA
  pAA = pA^2
  pD.aa = pD / (paa + r*pAa + r*pAA) #P(D|aa)
  pD.Aa = r*pD.aa
  pD.AA = r*pD.aa
  pX.D = c(pD.aa, pD.Aa, pD.AA)*c(paa, pAa, pAA)/pD
  pX.Dbar = (1-c(pD.aa, pD.Aa, pD.AA))*c(paa, pAa, pAA) / (1-pD)
  x=c(0, 1, 2)
  y=matrix(c(pX.D*ncase, pX.Dbar*nctrl), nrow = 3, ncol = 2, byrow = FALSE)
  y<-round(y)
  fit = glm(y~x, family=binomial(link = "logit")) #Null Deviance is the noncentral param.
  fit$null.deviance
}
