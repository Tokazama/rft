#' Estimates the smoothness/fwhm of an image matrix
#'
#'
#' @param mat-image matrix used for the fitted analysis
#' @param mask-mask image for the image matrix
#' @param psdf-pooled standard deviation obtained from fitted analysis
#' @return Outputs the estimated fwhm and covariance matrix that was used to estimate it
#' @examples
#'	
#'	mask<-getMask(imglist[[1]])
#'	mat <- imageListToMatrix(imglist, mask)
#'	var1<-c(1:nrow(mat))
#'  lmfit <- lm(mat~var1)
#'	res<-residuals(lmfit)
#'  res2<-colSums(res^2)
#'  rdf<-lmfit$df.residual
#'  S2<-res2/rdf
#'  psd<-sqrt(S2)
#'	fwhm<-est.smooth(mat,mask,psd)
#'
#' @export est.smooth
est.smooth<-function(Sres,mask,psd){
	voxels<-ncol(Sres)
	subs<-nrow(Sres)
	dimx<-dim(mask)[1]
	dimy<-dim(mask)[2]
	dimz<-dim(mask)[3]
	fwhm<-matrix(0L,nrow=1,ncol=3)
	
#Estimate partial derivatives of Znot at a voxel in x, y, and z direction
	partial.derivative<-function(img,x,y,z,lambda){	
		vox<-getPixels(img,x,y,z)[1]
		if(!is.na(vox)){
			xvox<-getPixels(img,x+1,y,z)[1]
			if(!is.na(xvox)){
				fwhm[1]<-fwhm[1]+(((xvox-vox)^2)/(voxels*(subs-1)))
				}
			yvox<-getPixels(img,x,y+1,z)[1]
			if(!is.na(yvox)){
				fwhm[2]<-fwhm[2]+(((yvox-vox)^2)/(voxels*(subs-1)))
				}
			zvox<-getPixels(img,x,y,z+1)[1]
			if(!is.na(zvox)){
				fwhm[3]<-fwhm[3]+(((zvox-vox)^2)/(voxels*(subs-1)))
				}
			}
		return(fwhm)
		}

#Estimate partial derivatives of standardized residuals from the fitted model
	fwhm2<-matrix(0L,nrow=1,ncol=3)
	for (i in 1:subs){
		img<-makeImage(mask,Sres[i,])
		for (x in 1:(dimx)){
			for (y in 1:(dimy)){
				for (z in 1:(dimz)){
					lambda<-partial.derivative(img,x,y,z,lambda)
					fwhm2<-fwhm2+fwhm
					}
				}
			}
		}


	vintegral<-function(t){((((t^2)+subs-1)^2)/((v-1)*(v-2)))*((dt(t,df)^3)/(dnorm(t,df)^2))}
	lamv<-integrate(vintegral,Inf,Inf)
	fwhm2<-sqrt(4*log(2)/fwhm2)*lamv
	return(fwhm2)
	}
