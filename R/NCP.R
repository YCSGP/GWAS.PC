#' This function calculates non-centrality parameter for single marker analysis.
#' @param r relative risk
#' @param pA minor allele frequency
#' @param ncase number of cases 
#' @param nctrl number of controls
#' @param pD disease prevalence
#' @param error.mat1 3*3 matrix of genotyping error in cases
#' @param error.mat2 3*3 matrix of genotyping error in controls
#' @param type "ndiff"/"diff", indicating non-differential or differential genotyping error
#' @param moi "m"/"r"/"d", indicating multiplicative, or recessive, or dominant mode of inheritance
#' @return the non-centrality parameter of the chisquare distribution 
#' @export
NCP.single<-function(r,pA, ncase, nctrl, pD, emat1, emat2=emat1, type="ndiff",moi="m")
{
  ncp<-numeric(1)
  if(type=="ndiff" && moi =="m" )  ncp<-mul.ndiff(r,pA, ncase, nctrl, pD, emat1)
  if(type=="ndiff" && moi =="r" )  ncp<-recessive.ndiff(r,pA, ncase, nctrl, pD, emat1)
  if(type=="ndiff" && moi =="d" )  ncp<-dominant.ndiff(r,pA, ncase, nctrl, pD, emat1)
  if(type=="diff"  && moi =="m")   ncp<-mul.diff(r,pA, ncase, nctrl, pD, emat1,emat2)
  if(type=="diff"  && moi =="r")   ncp<-recessive.diff(r,pA, ncase, nctrl, pD, emat1,emat2)
  if(type=="diff"  && moi =="d")   ncp<-dominant.diff(r,pA, ncase, nctrl, pD, emat1,emat2)
  ncp
}
#' This function calculates non-centrality parameter for multiple marker analysis.
#' @param r relative risk
#' @param pA minor allele frequency
#' @param ncase number of cases (assume number of controls is the same as cases)
#' @param pD disease prevalence
#' @param M number of markers in a gene/region
#' @param m10.case 1-sensitivity in cases
#' @param m00.case specificity in cases
#' @param m10.ctrl 1-sensitivity in controls
#' @param m00.ctrl specificity in controls
#' @param type "ndiff"/"diff", indicating non-differential or differential genotyping error
#' @return the non-centrality parameter of the chisquare distribution 
#' @export
NCP.multi<-function(r,pA, ncase, pD, M, m10.case, m00.case, m10.ctrl=m10.case, m00.ctrl=m00.case, type="ndiff")
{
  ncp<-numeric(1)
  if(type=="ndiff") ncp <- NCP(r,pA,pD,m10.case,m00.case,M,ncase)
  if(type=="diff")  ncp <- NCP.diff.alt(r,pA,pD,m10.case,m00.case,m10.ctrl,m00.ctrl,M,ncase)
  ncp
}


