#' This function returns the power of chi-square test given type I error and non-centrality parameter
#' @param alpha level of type I error
#' @param alt.ncp non-centrality parameter under alternative hypothesis
#' @param null.ncp non-centrality parameter under null hypohtesis
#' @return power of chi-square test
power.gwas<-function(alpha,alt.ncp,null.ncp=0)
{
  critical.t<-qchisq(1-alpha,df=1,ncp=null.ncp)
  gwas.power <- 1-pchisq(critical.t,df=1,ncp=alt.ncp)
  gwas.power
}
#' This function calculates the non-centrality parameter of the null distribution of single marker test
#' in existence of differential genotyping error
#' @param pA minor allele frequency of the SNP
#' @param ncase number of cases in the study
#' @param nctrl number of contrls in the study
#' @param error.mat.case 3*3 matrix of genotyping error in cases
#' @param error.mat.ctrl 3*3 matrix of genotyping error in controls
#' @return non-centrality parameter of the null distribution
diff.null<-function(pA, ncase, nctrl, error.mat.case,error.mat.ctrl)
{
  pa  = 1-pA 
  paa = pa^2 # P(aa)
  pAa = 2*pa*pA
  pAA = pA^2
  r = 1
  pX.D = c(pAA, pAa, paa)
  pX.Dbar = c(pAA, pAa, paa)
  pX.D=pX.D%*%error.mat.case
  pX.Dbar=pX.Dbar%*%error.mat.ctrl
  x=c(2, 1, 0)
  y=matrix(c(pX.D*ncase, pX.Dbar*nctrl), nrow = 3, ncol = 2, byrow = FALSE)
  y<-round(y)
  fit = glm(y~x, family=binomial(link = "logit")) #Null Deviance is the noncentral param.
  fit$null.deviance
}
#' This function calculates the statistical power for single marker analysis.
#' @param alpha type I error 
#' @param r relative risk
#' @param pA minor allele frequency
#' @param ncase number of cases 
#' @param nctrl number of controls
#' @param pD disease prevalence
#' @param error.mat1 3*3 matrix of genotyping error in cases
#' @param error.mat2 3*3 matrix of genotyping error in controls
#' @param type "ndiff"/"diff", indicating non-differential or differential genotyping error
#' @param moi "m"/"r"/"d", indicating multiplicative, or recessive, or dominant mode of inheritance
#' @return the statistical power of single marker association test
#' @export
Power.single<-function(alpha,r,pA, ncase, nctrl, pD, emat1, emat2=emat1, type="ndiff",moi="m")
{
  power1 <- numeric(1)
  ncp1<-NCP.single(r,pA, ncase, nctrl, pD, emat1, emat2=emat1, type,moi)
  if(type=="ndiff") power1 <- power.gwas(alpha,ncp1)
  if(type=="diff") {
    ncp0 <- diff.null(pA, ncase, nctrl, emat1,emat2)
    power1 <- power.gwas(alpha,alt.ncp=ncp1,null.ncp=ncp0)
  }
  power1
}
#' This function calculates the statistical power for multiple marker analysis.
#' @param alpha type I error 
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
Power.multi<-function(alpha,r,pA, ncase, pD, M, m10.case, m00.case, m10.ctrl=m10.case, m00.ctrl=m00.case, type="ndiff")
{
  power1<-numeric(1)
  ncp1<-NCP.multi(r,pA, ncase, pD, M, m10.case, m00.case, m10.ctrl, m00.ctrl, type)
  
  if(type=="ndiff") power1 <- power.gwas(alpha,ncp1)
  if(type=="diff") {
    ncp0 <- NCP.diff.null(pA,pD,m10.case,m00.case,m10.ctrl,m00.ctrl,M,ncase)
    power1 <- power.gwas(alpha,alt.ncp=ncp1,null.ncp=ncp0)
  }
  power1
}


