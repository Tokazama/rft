#' RFT p-values
#'
#' This function calculates the probability of obtaining a cluster the size of 
#' code/{cMask} within a search region of code/{bMask} given the statistical threshold
#' of code/{stat} used to extract said cluster. The degrees of freedom and statistical 
#' field type used to obtain the original statistical map are also required. 
#' 
#' set-level
#' rft.pval(D, c, k, u, resels, df, fieldType)
#' 
#' cluster-level
#' rft.pval(D, 1, k, u, resels, df, fieldType)
#' 
#' voxel-level
#' rft.pval(D, 1, 0, u, resels, df, fieldType)
#' 
#' @param D image dimensions
#' @param c Threshold
#' @param k Cluster size in resel space
#' @param u Number of clusters
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @param fieldType:
#' \describe{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param n number of statistical fields in conjunction
#' 
#' @return The probability of obtaining the specified cluster
#' Pcor corrected p-value
#' Pu uncorrected p-value
#' Ec expected number of clusters
#' ek expected number of resels per cluster
#' @reference 
#' Friston K.J., (1994) Assessing the Significance of Focal Activations Using Their Spatial Extent.
#' Friston K.J., (1996) Detecting Activations in PET and fMRI: Levels of Inference and Power.
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' @Author Zachary P. Christensen
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
#' @examples
#' 
#' 
#' 
#' @export rft.pval
rft.pval <-function(D, c, k, u, resels, df, fieldType, n=1){
  if(missing(fieldType)){
    stop("Must specify fieldType")
  }else if(missing(df)){
    stop("Must specify df")
  }else if(missing(resels)){
    stop("Must specify resels")
  }else if(missing(u) && missing(k) && missing(c)){
    stop("Must atleast specify one of u, k, or c")
  }
  G <-sqrt(pi)/gamma((1:(D+1)/2))
  ec <-c(0,0,0,0)
  if (fieldType=="T"){
    ec[1] <-1-pt(u,df[2])
    ec[2] <-(((4*log(2))^(1/2))/(2*pi))*((1+((u^2)/df[2]))^(-1/2*(df[2]-1)))
    ec[3] <-(4*log(2))/((2*pi)^(3/2))*((1+u^2/df[2])^((1-df[2])/2))*u/((df[2]/2)^(1/2))*exp(lgamma((df[2]+1)/2)-lgamma(df[2]/2))
    ec[4] <-(((4*log(2))^(3/2))/((2*pi)^2))*((1+((u^2)/df[2]))^(-1/2*(df[2]-1)))*((((df[2]-1)/df[2])*(u^2))-1)
  }else if(fieldType=="F"){
    ec[1] <-1-pf(u,df[1],df[2])
    ec[2] <-((4*log(2))/(2*pi))^(1/2)*exp(lgamma((df[2]+df[1]-1)/2)-(lgamma(df[2]/2) + lgamma(df[1]/2)))*2^(1/2)*(df[1]*u/df[2])^(1/2*(df[1]-1))*(1+df[1]*u/df[2])^(-1/2*(df[2]+df[1]-2))
    ec[3] <-((4*log(2))/(2*pi))*exp(lgamma((df[2]+df[1]-2)/2)-(lgamma(df[2]/2) + lgamma(df[1]/2)))*(df[1]*u/df[2])^(1/2*(df[1]-2))*(1+df[1]*u/df[2])^(-1/2*(df[2]+df[1]-2))*((df[2]-1)*df[1]*u/df[2]-(df[1]-1))
    ec[4] <-((4*log(2))/(2*pi))^(3/2)*exp(lgamma((df[2]+df[1]-3)/2)-(lgamma(df[2]/2) + lgamma(df[1]/2)))*2^(-1/2)*(df[1]*u/df[2])^(1/2*(df[1]-3))*(1+df[1]*u/df[2])^(-1/2*(df[2]+df[1]-2))*((df[2]-1)*(df[2]-2)*(df[1]*u/df[2])^2-(2*df[2]*df[1]-df[2]-df[1]-1)*(df[1]*u/df[2])+(df[1]-1)*(df[1]-2))
  }else if(fieldType=="X"){
    ec[1] <-1-pchisq(u,df[2])
    ec[2] <-((4*log(2))/(2*pi))^(1/2)*(u^(1/2*(df[2] - 1))*exp(-u/2-lgamma(df[2]/2))/2^((df[2]-2)/2))
    ec[3] <-((4*log(2))/(2*pi))*(u^(1/2*(df[2] - 1))*exp(-u/2-lgamma(df[2]/2))/2^((df[2]-2)/2))*(u-(df[2]-1))
    ec[4] <-((4*log(2))/(2*pi))^(3/2)*(u^(1/2*(df[2] - 1))*exp(-u/2-lgamma(df[2]/2))/2^((df[2]-2)/2))*(u^2-(2*df[2]-1)*u+(df[2]-1)*(df[2]-2))
  }else if(fieldType=="Z"){
    ec[1] <-1-pnorm(u,df[2])
    ec[2] <-(4*log(2))^(1/2)/(2*pi)*exp(-u^2/2)
    ec[3] <-(4*log(2))/((2*pi)^(3/2))*exp(-u^2/2)*u
    ec[4] <-(4*log(2))^(3/2)/((2*pi)^2)*exp(-u^2/2)*(u^2 - 1)
  }
  ec <-pmax(ec[1:(D+1)],.Machine$double.eps)
  P <-toeplitz(as.numeric(ec*G))
  P[lower.tri(P)] <-0
  if (n != round(n)){
    n <- round(n)
    warning("rounding exponent `n' to", n)
  }
  phi <-diag(nrow=nrow(P))
  pot <-P
  while (n > 0) {
    if (n %% 2) 
      phi <-phi %*% pot
    n <-n%/%2
    pot <-pot %*% pot
  }
  P <-phi
  P <-P[1,]
  EM <-(resels[1:(D+1)]/G)*P # maxima in all dimensions
  Ec <-sum(EM) # number of overall expected maxima/clusters
  EN <-P[1]*resels[D+1] # number of resels in entire image
  ek <-EN/EM[D+1] # expected number of resels per cluster
  
  rfB <-(gamma(D/2+1)/ek)^(2/D)
  Pu <-exp(-rfB*(k^(2/D))) # cumulative cluster-size distribution from which uncorrected P values are calculated
  
  Pcor <-1-ppois(c-1,lambda=(Ec+.Machine$double.eps)*Pu) 
  z <-list(Pcor=Pcor, Pu=Pu, Ec=Ec, ek=ek,ec=ec)
  return(z)
}
