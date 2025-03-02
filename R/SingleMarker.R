#' This function is called in burden test
#' @param r relative risk
#' @param pA allele frequency
#' @param pD disease prevalence
#' @return genotype distribution break down by cases and controls
SingleMarker<-function(r,pA,pD){
  pa  = 1-pA  ### Attention a is Major and A is minor
  paa = pa^2 # P(aa)
  pAa = 2*pa*pA
  pAA = pA^2
  pD.aa = pD / (paa + r*pAa + r^2*pAA) #P(D|aa)
  pD.Aa = r*pD.aa
  pD.AA = r^2*pD.aa
  pX.D = c(pD.aa, pD.Aa, pD.AA)*c(paa, pAa, pAA)/pD
  pX.Dbar = (1-c(pD.aa, pD.Aa, pD.AA))*c(paa, pAa, pAA) / (1-pD)
  p01<-pX.D[1]
  p11<-pX.D[2]+pX.D[3]
  p00<-pX.Dbar[1]
  p10<-pX.Dbar[2]+pX.Dbar[3]
  c(p01,p11,p00,p10)
}
