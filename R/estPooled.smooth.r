estPooled.smooth<-function(res,rdf,mask){
	res2<-colSums(res^2)
	rdf<-regfit$df.residual
	S2<-res2/rdf
	psd<-sqrt(S2)
	Mmat<-colMeans(res)
	Zmat<-res
	fwhm<-matrix(0L,nrow=1,ncol=3)
	subs<-nrow(res)

	cat("Estimating fwhm/smoothing")
	progress <- txtProgressBar(min = 0, max = subs, style = 3)
	for (i in 1:subs){
		Zmat[i,]<-(res[i,]-Mmat[1])/psd
		img<-makeImage(mask,Zmat[i,])
		smooth<-est.smooth(img,mask,1,1,1)
		fwhm<-fwhm+smooth[[2]]
		setTxtProgressBar(progress, i)
		}
	close(progress)
	fwhm2<-sqrt(4*log(2)/(fwhm/(subs-1)))
	return(fwhm2)
	}
