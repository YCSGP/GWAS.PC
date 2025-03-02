#' This function calculates the non-centrality parameter of burden test statistic under null hypothesis
#' when the genotyping error is different in cases and controls.
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param m10.case 1-sensitivity in cases
#' @param m00.case specificity in cases
#' @param m10.ctrl 1-sensitivity in controls
#' @param m00.ctrl specificity in controls
#' @param M number of markers in a gene/region
#' @param Ncase number of cases 
#' @return the non-centrality parameter of burden test statistic under null hypothesis
#' @export
NCP.diff.null<-function(pA,pD,m10.case,m00.case,m10.ctrl,m00.ctrl,M,Ncase)
{
  phi.case<-PhiCaseWES(1,pA,pD,m10.case,m00.case,M)
  phi.ctrl<-PhiCtrlWES(pA,pD,m10.ctrl,m00.ctrl,M)
  ncp<-Ncase*((phi.case-phi.ctrl)^2/(phi.case+phi.ctrl)+(phi.case-phi.ctrl)^2/(2-phi.case-phi.ctrl))
  ncp
}
#' This function calculates the non-centrality parameter of burden test statistic under alternative hypothesis
#' when the genotyping error is different in cases and controls.
#' @param r relative risk
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param m10.case 1-sensitivity in cases
#' @param m00.case specificity in cases
#' @param m10.ctrl 1-sensitivity in controls
#' @param m00.ctrl specificity in controls
#' @param M number of markers in a gene/region
#' @param Ncase number of cases 
#' @return the non-centrality parameter of burden test statistic under alternative hypothesis
#' @export
NCP.diff.alt<-function(r,pA,pD,m10.case,m00.case,m10.ctrl,m00.ctrl,M,Ncase)
{
  phi.case<-PhiCaseWES(r,pA,pD,m10.case,m00.case,M)
  phi.ctrl<-PhiCtrlWES(pA,pD,m10.ctrl,m00.ctrl,M)
  ncp<-Ncase*((phi.case-phi.ctrl)^2/(phi.case+phi.ctrl)+(phi.case-phi.ctrl)^2/(2-phi.case-phi.ctrl))
  ncp
}
