thresh.rft<-function(FWHM,k,D,df,fieldmat,fieldtype=T,mask){
  #S-voxels in image
  #k-expected voxels in uperthreshold of cluster
  #alpha-p value
  #
  S<-labelGeometryMeasures(mask)[2,2]
  u<-1
  EN<-S*alpha
  if (EN < k){
    stop("Expected cluster size is too small")
  }
  ka<-k-1
  if (fieldtype=="T"){
    while(ka > k){
      EN<-S*(1-pt(u,df))
      W<-FWHM/((4*log(2))^(1/2))
      Em<-S*(2*pi)^(-(D+1)/2)*(W^(-D))*(u^(D-1))*(exp((-u^2)/2))
      beta<-(gamma((D/2)+1)*Em/EN)^(2/D)
      ka<-(log(-Em/log(1-alpha))/beta)^(D/2)
      u<-u+.01
    }
  }else if(fieldtyp="F"){
    while(ka > k){
      EN<-S*(1-pf(u, df[1],df[2]))
      W<-FWHM/((4*log(2))^(1/2))
      Em<-S*(2*pi)^(-(D+1)/2)*(W^(-D))*(u^(D-1))*(exp((-u^2)/2))
      beta<-(gamma((D/2)+1)*Em/EN)^(2/D)
      ka<-(log(-Em/log(1-alpha))/beta)^(D/2)
      u<-u+.01
    }
  }else if(fieldtyp="X"){
    while(ka > k){
      EN<-S*(1-pchisq(u, df[1],df[2]))
      W<-FWHM/((4*log(2))^(1/2))
      Em<-S*(2*pi)^(-(D+1)/2)*(W^(-D))*(u^(D-1))*(exp((-u^2)/2))
      beta<-(gamma((D/2)+1)*Em/EN)^(2/D)
      ka<-(log(-Em/log(1-alpha))/beta)^(D/2)
      u<-u+.01
    }
  }else if(fieldtyp="G"){
    while(ka > k){
      EN<-S*(1-qnorm(u))
      W<-FWHM/((4*log(2))^(1/2))
      Em<-S*(2*pi)^(-(D+1)/2)*(W^(-D))*(u^(D-1))*(exp((-u^2)/2))
      beta<-(gamma((D/2)+1)*Em/EN)^(2/D)
      ka<-(log(-Em/log(1-alpha))/beta)^(D/2)
      u<-u+.01
    }
  }else {
    stop("Error in calculating threshold for image")
  }
  img<-makeImage(mask,fieldmat)
  thrimg<-thresholdImage(img,u,Inf)
  clusters<-labelClusters(thrimg)
}
