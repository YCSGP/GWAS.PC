#' This function calculates the proportion of subjects with at least one rare allele in cases
#' when the genotyping error is the same in cases and controls.
#' @param r relative risk
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param m10 1-sensitivity
#' @param m00 specificity
#' @param M number of markers in a gene/region
#' @return the proportion of subjects with at least one rare allele in cases
PhiCaseWES<-function(r,pA,pD,m10,m00,M)
{
  sm<-SingleMarker(r,pA,pD/M)
  p01<-sm[1];p11<-sm[2];p00<-sm[3];p10<-sm[4];
  sum.vec<-numeric(M+1)
  Pr1<-numeric(M+1)
  for(k in 0:M) Pr1[k+1]<-k/M*p11*p10^(k-1)*p00^(M-k)+(M-k)/M*p01*p00^(M-k-1)*p10^k
  for(k in 0:M) sum.vec[k+1]<-choose(M,k)*m00^(M-k)*m10^k*Pr1[k+1]
  phi.case.wes<-1-sum(sum.vec)
  phi.case.wes
}
#' This function calculates the proportion of subjects with at least one rare allele in controls
#' when the genotyping error is the same in cases and controls.
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param m10 1-sensitivity
#' @param m00 specificity
#' @param M number of markers in a gene/region
#' @return the proportion of subjects with at least one rare allele in cases
PhiCtrlWES<-function(r,pA,pD,m10,m00,M)
{
  sm<-SingleMarker(r,pA,pD/M)
  p01<-sm[1];p11<-sm[2];p00<-sm[3];p10<-sm[4];
  sum.vec<-Pr1<-numeric(M+1)
  for(k in 0:M) Pr1[k+1]<-p10^k*p00^(M-k)
  for(k in 0:M) sum.vec[k+1]<-choose(M,k)*m10^k*m00^(M-k)*Pr1[k+1]
  phi.ctrl.wes<-1-sum(sum.vec)
  phi.ctrl.wes
}
#' This function calculates the non-centrality parameter of burden test statistic under alternative hypothesis
#' when the genotyping error is the same in cases and controls.
#' @param r relative risk
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param m10 1-sensitivity
#' @param m00 specificity
#' @param M number of markers in a gene/region
#' @param Ncase number of cases 
#' @return the non-centrality parameter of burden test statistic under alternative hypothesis
#' @export
NCP<-function(r,pA,pD,m10,m00,M,Ncase)
{
  phi.case<-PhiCaseWES(r,pA,pD,m10,m00,M)
  phi.ctrl<-PhiCtrlWES(r,pA,pD,m10,m00,M)
  ncp<-Ncase*((phi.case-phi.ctrl)^2/(phi.case+phi.ctrl)+(phi.case-phi.ctrl)^2/(2-phi.case-phi.ctrl))
  ncp
}

