ants.rft<-function(D,k,lmfit,mat,fieldtype,mask){
	rdf<-lmfit$df.residual
	res<-residuals(lmfit)
	sumres<-(rowSums(res)^2)/rdf
	sqres<-rowSums(res^2)
	ssq<-sumres-sqres
	for (i in 1:nrow(res)){
		res[i,]<-res[i,]/sqrt(ssq[i,])
		}
	fwhm<-est.smooth(Sres,mask,psdf)
	thresh.imgs<-thresh.rft(D,fwhm,k,df,stat,fieldtype,mask)
	limg<-length(thresh.img)
	for (i in 1:limg){
		res<-paste("res",i, sep="")
    res<-ants.resels(thesh.img[i],fwhm)
    cluster.fwhm<-est.smooth(mat,tresh.img[i],psd)
    pval[i]<-ants.ec(mat,fieldtype,df,res)
  }
  
   return(pval
}
