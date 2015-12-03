#' @name rft.ec
#' @title Calculates the euler characteristic
#'
#'
#' @param stat - statistical value (typically the maxima of a cluster or SPM)
#' @param fieldType:
#'	"T"- T-field
#'	"F"- F-field
#'	"X"- Chi squar field
#'	"Z"- Gaussian field
#' @param df - degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @return ec - a vector of estimated euler characteristics
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
rft.ec<-function(stat,fieldType,df){
  ec<-c(0,0,0,0)
  t<-stat
  if(fieldType=="T"){
    a<-gamma((df+1)/2)/(((df*pi)^(1/2))*gamma(df/2))
    tintegral<-function(stat){a*((1+(stat^2)/df)^(-1/2*(df+1)))}
    iec<-integrate(tintegral,lower=stat,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/(2*pi))*((1+((stat^2)/df))^(-1/2*(df-1)))
    ec[3]<-((((4*log(2))*lgamma((df+1)))/(((2*pi)^(3/2))*((df/2)^(1/2))*lgamma(df/2)))*((1+((stat^2)/2))^(-1/2*(df-1))))*stat
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^2))*((1+((stat^2)/df))^(-1/2*(df-1)))*((((df-1)/df)*(stat^2))-1)
  }else if(fieldType=="F"){
    a<-(gamma((df+k-2)/2)/(gamma(df/2)*gamma(k/2)))
    fintegral<-function(stat){a*(((k*stat/df)^(1/2*(k-2)))*(1+(k*stat/2)^(-1/2*(df+k))))}
    iec<-integrate(fintegral,lower=stat,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/((2*pi)^(1/2)))*(((lgamma(df+k-1)/2)*2^(1/2))/(lgamma(df/2)*lgamma(k/2)))*((k*stat/df)^(1/2(k-2)))*((1+(k*stat/v))^(-1/2*(df+k-2)))
    ec[3]<-((4*log(2))/(2*pi))*((lgamma((v+k-2)))/(lgamma(df/2)*lgamma(k/2)))*((k*stat/df)^(1/2*(k-2)))*(1+(k*t/df))^(-1/2(df+k-2))*((v-1)*(k*stat/v)-(k-1))
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^(3/2)))*(((lgamma((v+k-2)))*2^(-1/2))/(lgamma(df/2)*lgamma(k/2)))*((k*t/df)^(1/2*(k-3)))*((1+(k*stat/df))^(-1/2*(df+k-2)))*((df-1)*(df-2)*((k*stat/v)^(2))-(2*df*k-df-k-1)*(k*stat/df)+(k-1)*(k-2))
  }else if(fieldType=="X"){
    xintegral<-function(stat){((stat^(1/2))*e^(-1/2*stat))/((2^(df/2))*(gamma(df/2)))}
    iec<-integrate(xintegral,lower=stat,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-(((4*log(2))^(1/2))/(2*pi))*(((stat^(1/2*(df-1)))*(e^(-1/2*stat)))/((2^(1/2(df-2)))*gamma(df/2)))
    ec[3]<-((4*log(2))/(2*pi))*(((stat^(1/2*(df-2)))*(e^(-1/2*stat)))/((2^(1/2*(df-2)))*gamma(df/2)))*(stat-(v-1))
    ec[4]<-((4*log(2)^(3/2))/((2*pi)^(3/2)))*(((stat^(1/2*(df-3)))*(e^(-1/2*stat)))/((2^(1/2*(df-2)))*(gamma(df/2))))*((stat^2)-((2*df-1)*stat)+(df-1)*(df-2))
  }else if(fieldType=="G"){  
    gintegral<-function(stat){((1/((2*pi)^(1/2)))*e^(-1/2*(stat^2)))}
    iec<-integrate(gintegral,lower=stat,upper=Inf)
    ec[1]<-as.numeric(iec[1])
    ec[2]<-((((4*log(2))^(1/2))/(2*pi))*e^(-1/2*((stat^2))))
    ec[3]<-((4*log(2))/(2*pi))*stat*(e^(-1/2*(stat^2)))
    ec[4]<-(((4*log(2))^(3/2))/((2*pi)^2))*((stat^2)-1)*(e^(-1/2*(stat^2)))
  }else{
    stop("Incorrect input for FieldType")
    }
  return(ec)
}
