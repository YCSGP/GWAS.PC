#' This function calculates the minimum sample size needed to achieve a specified level 
#' of power for single marker analysis.
#' @param alpha type I error 
#' @param r relative risk
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param error.mat1 3*3 matrix of genotyping error in cases
#' @param error.mat2 3*3 matrix of genotyping error in controls
#' @param type "ndiff"/"diff", indicating non-differential or differential genotyping error
#' @param moi "m"/"r"/"d", indicating multiplicative, or recessive, or dominant mode of inheritance
#' @param p pre-specified statistical power
#' @return the minimum sample size needed to achieve a specified level 
#' @return of power for single marker analysis.
#' @export
Size.single<-function(alpha,p, r,pA, pD, emat1, emat2=emat1, type="ndiff",moi="m",range=c(1,100000))
{
  ivec<-seq(range[1],range[2],by=10)
  power1<-numeric(length(ivec))
  for(i in 1:length(ivec))
    power1[i] <- Power.single(alpha,r,pA, ivec[i], ivec[i], pD, emat1, emat2, type,moi)
  if(max(power1) <= p) {print("Need to expand range!"); return(0); } else {
    imin<-min(which(power1>=p))
    return(ivec[imin])
  }
}
#' This function calculates the minimum sample size needed to achieve a specified level 
#' of power for multiple marker analysis.
#' @param alpha type I error 
#' @param r relative risk
#' @param pA minor allele frequency
#' @param pD disease prevalence
#' @param M number of markers in a gene/region
#' @param m10.case 1-sensitivity in cases
#' @param m00.case specificity in cases
#' @param m10.ctrl 1-sensitivity in controls
#' @param m00.ctrl specificity in controls
#' @param type "ndiff"/"diff", indicating non-differential or differential genotyping error
#' @param p pre-specified statistical power
#' @return the minimum sample size needed to achieve a specified level 
#' @return of power for multiple marker analysis.
#' @export
Size.multi<-function(alpha,p,r,pA, pD, M, m10.case, m00.case, m10.ctrl=m10.case, 
                      m00.ctrl=m00.case, type="ndiff",range=c(1,100000))
{
    ivec<-seq(from=range[1],to=range[2],by=10)
    power1<-numeric(length(ivec))
    for(i in 1:length(ivec))
        power1[i] <- Power.multi(alpha,r,pA, ivec[i], pD, M, m10.case, m00.case, m10.ctrl, m00.ctrl, type)
    if(max(power1) <= p) {print("Need to expand range!"); return(0); } else {
    imin<-min(which(power1>=p))}
    return(ivec[imin])
}
  

