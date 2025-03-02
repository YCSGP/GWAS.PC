#' This function returns the non-centrality parameter when the genotyping error is different in cases
#' and in controls, assuming multiplicative mode of inheritance.
#' @param r relative risk
#' @param pA minor allele frequency
#' @param ncase number of cases
#' @param nctrl number of control samples
#' @param pD disease prevalence
#' @param error.mat.case 3*3 matrix of genotyping error in cases
#' @param error.mat.ctrl 3*3 matrix of genotyping error in controls
#' @return non-centrality parameter of the chi-square distribution under alternative hypothesis
#' @export
mul.diff<-function(r,pA, ncase, nctrl, pD, error.mat.case,error.mat.ctrl)
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
  pX.D=pX.D%*%error.mat.case
  pX.Dbar=pX.Dbar%*%error.mat.ctrl
  x=c(2, 1, 0)
  y=matrix(c(pX.D*ncase, pX.Dbar*nctrl), nrow = 3, ncol = 2, byrow = FALSE)
  y<-round(y)
  fit = glm(y~x, family=binomial(link = "logit")) #Null Deviance is the noncentral param.
  fit$null.deviance
}
