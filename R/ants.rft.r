ants.rft<-function(D,k,lmfit,mat,fieldtype,mask){
	res<-residuals(lmfit)
	res2<-colSums(res^2)
	rdf<-lmfit$df.residual
	S2<-res2/rdf
	psd<-sqrt(S2)
	
	
	fwhm<-est.smooth(mat,mask,psdf)
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
