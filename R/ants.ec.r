#' @name rft.ec
#' @title Calculates the euler character 
#'
#'
#' @param stat-statistical value (typically the maxima of a cluster or SPM)
#' @param fieldType:
#'	"T"- T-field
#'	"F"- F-field
#'	"X"- Chi squar field
#'	"Z"- Gaussian field
#' @param df- degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @return ec-a vector of estimated euler characteristics
#'	
#'  var1<-vardata[,10]
#'  subs<-nrow(varmat)
#'  voxels<-ncol(varmat)
#'  regpval<-matrix(nrow=1,ncol=voxels)
#'  regtstat<-matrix(nrow=1,ncol=voxels)
#'  resmat<-matrix(OL,nrow=subs,ncol=voxels)
#'  for (i in 1:voxels){
#'	  vox<-varmat[,i]
#'	  regfit<-lm(vox~var1)
#'	  ##Extract statistical values
#'	  resmat[,i]<-residuals(regfit)
#'	  regsum<-summary(regfit)
#'	  regtstat[,i]<-regsum$coefficients[3,3]
#'	  }
#'  fwhm<-estPooled.smooth(res,rdf,mask)
#'	negclust<-image2ClusterImages(timg,150,-Inf,-thresh)
#'	negtable<-matrix(ncol=5)
#'	colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
#'	for (i in 1:length(posclust)){
#'	  cat("Determing negative cluster level statistics:",i,sep=" ")
#'	  cmask<-getMask(negclust[[i]])
#'	  cvoxs<-sum(as.array(cmask))
#'	  pclust<-rft.pcluster(negclust[[i]],mask,fwhm,thresh,df,fieldType)
#'	  peak<-max(posclust[[i]])
#'	  loc<-labelImageCentroids(posclust[[i]])[2]
#'	  resel <-ants.resel(mask,fwhm)
#'	  ec<-ants.ec(stat,fieldType,df)
#'	  pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
#'	  negtable[i,]<-c(cvox, pclust, peak, pval, loc$vertices[1],loc$vertices[2],loc$vertices[3])
#'	  clustername<-paste("N-Cluster:",i,sep="")
#'	  rownames(postable[i,])<-c(clustername)
#'	  image<-paste(fileDir,"Nlcuster",i,".nii.gz",sep="")
#'	  antsImageWrite(negclust[[i]],file=image)
#'	  }
#'
#' @export ants.ec
ants.ec<-function(stat,fieldType,df){
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
