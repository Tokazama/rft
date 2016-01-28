#' RFT p-values
#'
#' This function calculates the probability of obtaining a cluster the size of 
#' code/{cMask} within a search region of code/{bMask} given the statistical threshold
#' of code/{stat} used to extract said cluster. The degrees of freedom and statistical 
#' field type used to obtain the original statistical map are also required. 
#' 
#' set-level
#' rft.pval(c, k, u, resels, df, fieldType)
#' 
#' cluster-level
#' rft.pval(1, k, u, resels, df, fieldType)
#' 
#' voxel-level
#' rft.pval(1, 1, u, resels, df, fieldType)
#' 
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
#' ## estimatation of a single images smoothness
#' outimg1 <-makeImage(c(10,10,10), rnorm(1000))
#' maskimg <-getMask(outimg1)
#' fwhm <-est.Smooth(outimg1,maskimg)
#' ## create arbitrary threshold and degrees of freedom from hypothetical analysis
#' thresh <-mean(outimg1)
#' df <-4
#' clustimg <-thresholdImage(outimg1, thresh, Inf, inval=1, outval=0)
#' pval <-rft.pcluster(3, clustmask, clustimg, thresh, df, fwhm[[1]], fieldType="T")
#'
#' @export rft.pval
rft.pval <-function(D, c, k, u, resels, df, fieldType){
  if(missing(fieldType)){
    stop("Must specify fieldType")
  }else if(missing(df)){
    stop("Must specify df")
  }else if(missing(resels)){
    stop("Must specify resels")
  }else if(missing(u) && missing(k) && missing(c)){
    stop("Must atleast specify one of u, k, or c")
  }

  n <-1
  
  G <-sqrt(pi)/gamma((1:(D+1)/2))
  euler <-rft.ec(u, fieldType, df[2])
  euler <-pmax(euler[1:(D+1)],.Machine$double.eps)
  P <-toeplitz(as.numeric(euler*G)) 
  P[lower.tri(P)] <-0
  # P <-mx.exp(P, n)
  P <-P[1,]
  EM <-(resels[1:(D+1)]/G)*P # maxima in all dimensions
  Ec <-sum(EM) # number of overall expected maxima/clusters
  EN <-P[1]*resels[D+1] # number of resels
  ek <-EN/EM[D+1] # expected number of resels per cluster
  
  rfB <-(gamma((D/2)+1)*ek)^(2/D)
  Pu <-exp(-rfB*(k^(2/D)))
  
  if (c > 1){
    # set-level
    Pcor <-1-ppois(c-1,lambda=(Ec*Pu))
  }else if (c==1 && k > 1){
    # cluster-level
    Pcor <-1-exp(-Ec*Pu)
  }else if (c==1 && k==1){
    # voxel-level
    Pcor <-1-exp(-Ec)
  }
  z <-list(Pcor=Pcor, Pu=Pu, Ec=Ec, ek=ek)
  z
}

#' Calculates the euler characteristic
#'
#' Used in conjuction with \code{rft.resel} the voxel-level probability of a
#' single voxel value can be determined. The voxel value is typically the peak 
#' voxel of a cluster after thresholding a statistical image. 
#' 
#' 
#' @param u statistical value (typically the maxima of a cluster or statistical map)
#' @param fieldType:
#' \describe{
#' \item{"T"}{T-field} 
#' \item{"F"}{F-field} 
#' \item{"X"}{Chi-square field"} 
#' \item{"Z"}{Gaussian field}
#' }
#' @param df degrees of freedom expressed as df[degrees of interest, degrees of error]
#'
#' @return ec a vector of estimated euler characteristics
#'
#' @References Worlsey K.J., et al. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' @Author Zachary P. Christensen
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
#' @examples
#' 
#' outimg1 <-makeImage(c(10,10,10), rnorm(1000))
#' maskimg <-getMask(outimg1)
#' fwhm <-est.Smooth(outimg1, maskimg)
#' resel <-rft.resels(maskimg, fwhm[[1]])
#' ec <-rft.ec(max(outimg1), fieldType="T", rnorm(1))
#' pvox<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
#' 
#' 
#' @export rft.ec
rft.ec<-function(u,fieldType,df){
  ec<-c(0,0,0,0)
  t<-u
  if(fieldType=="T"){
    a<-gamma((df[2]+1)/2)/(((df[2]*pi)^(1/2))*gamma(df[2]/2))
    tintegral<-function(u){a*((1+(u^2)/df[2])^(-1/2*(df[2]+1)))}
    iec<-integrate(tintegral,lower=u,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/(2*pi))*((1+((u^2)/df[2]))^(-1/2*(df[2]-1)))
    ec[3]<-((((4*log(2))*lgamma((df[2]+1)))/(((2*pi)^(3/2))*((df[2]/2)^(1/2))*lgamma(df[2]/2)))*((1+((u^2)/2))^(-1/2*(df[2]-1))))*u
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^2))*((1+((u^2)/df[2]))^(-1/2*(df[2]-1)))*((((df[2]-1)/df[2])*(u^2))-1)
  }else if(fieldType=="F"){
    a<-(gamma((df+k-2)/2)/(gamma(df/2)*gamma(k/2)))
    fintegral<-function(u){a*(((k*u/df)^(1/2*(k-2)))*(1+(k*u/2)^(-1/2*(df+k))))}
    iec<-integrate(fintegral,lower=u,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/((2*pi)^(1/2)))*(((lgamma(df+k-1)/2)*2^(1/2))/(lgamma(df/2)*lgamma(k/2)))*((k*u/df)^(1/2(k-2)))*((1+(k*u/v))^(-1/2*(df+k-2)))
    ec[3]<-((4*log(2))/(2*pi))*((lgamma((v+k-2)))/(lgamma(df/2)*lgamma(k/2)))*((k*u/df)^(1/2*(k-2)))*(1+(k*t/df))^(-1/2(df+k-2))*((v-1)*(k*u/v)-(k-1))
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^(3/2)))*(((lgamma((v+k-2)))*2^(-1/2))/(lgamma(df/2)*lgamma(k/2)))*((k*t/df)^(1/2*(k-3)))*((1+(k*u/df))^(-1/2*(df+k-2)))*((df-1)*(df-2)*((k*u/v)^(2))-(2*df*k-df-k-1)*(k*u/df)+(k-1)*(k-2))
  }else if(fieldType=="X"){
    xintegral<-function(u){((u^(1/2))*e^(-1/2*u))/((2^(df/2))*(gamma(df/2)))}
    iec<-integrate(xintegral,lower=u,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/(2*pi))*(((u^(1/2*(df-1)))*(e^(-1/2*u)))/((2^(1/2(df-2)))*gamma(df/2)))
    ec[3]<-((4*log(2))/(2*pi))*(((u^(1/2*(df-2)))*(e^(-1/2*u)))/((2^(1/2*(df-2)))*gamma(df/2)))*(u-(v-1))
    ec[4]<-((4*log(2)^(3/2))/((2*pi)^(3/2)))*(((u^(1/2*(df-3)))*(e^(-1/2*u)))/((2^(1/2*(df-2)))*(gamma(df/2))))*((u^2)-((2*df-1)*u)+(df-1)*(df-2))
  }else if(fieldType=="G"){  
    gintegral<-function(u){((1/((2*pi)^(1/2)))*e^(-1/2*(u^2)))}
    iec<-integrate(gintegral,lower=u,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-((((4*log(2))^(1/2))/(2*pi))*e^(-1/2*((u^2))))
    ec[3]<-((4*log(2))/(2*pi))*u*(e^(-1/2*(u^2)))
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^2))*((u^2)-1)*(e^(-1/2*(u^2)))
  }else{
    stop("Incorrect input for FieldType")
  }
  return(ec)
}
