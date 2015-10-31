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
est.smooth<-function(mat,mask,psd){
	voxels<-ncol(mat)
	subs<-length(mat)
	Mmat<-colMeans(mat)
	dimx<-dim(mask)[1]
	dimy<-dim(mask)[2]
	dimz<-dim(mask)[3]
	for (i in 1:subs){
		Zmat[i,]<-(mat[i,]-Mmat)/psd
		}
	lambda<-matrix(0L,nrow=3,ncol=3)
	for (i in 1:subs){
		img<-makeImage(mask,Zmat[i,]
		for (x in 1:(dimx)){
			for (y in 1:(dimy)){
				for (z in 1:(dimz)){
					vox<-getPixels(img,x,y,z)[1]
						if(!is.na(vox){
							xvox<-getPixels(img,x+1,y,z)[1]
							yvox<-getPixels(img,x,y+1,z)[1]
							zvox<-getPixels(img,x,y,z+1)[1]
							if(!is.na(xvox){
								lambda[1,1]<-lambda[1,1]+(((xvox-vox)^2)/(voxels*(subs-1)))
								}
							if(!is.na(yvox){
								lambda[2,2]<-lambda[2,2]+(((yvox-vox)^2)/(voxels*(subs-1)))
								}
							if(!is.na(zvox){
								lambda[3,3]<-lambda[3,3]+(((-vox)^2)/(voxels*(subs-1)))
								}
							if(!is.na(xvox) && !is.na(yvox)){	
								lambda[1,2]<-lambda[1,2]+(((vox+getPixels(img,x,y+1,z)[1])*(vox+getPixels(img,x+1,y,z)[1]))/(4*voxels*(subs-1)))
								}
							if(!is.na(xvox) && !is.na(zvox){
								lambda[1,3]<-lambda[1,3]+(((vox+getPixels(img,x,y,z+1)[1])*(vox+getPixels(img,x+1,y,z)[1]))/(4*voxels*(subs-1)))
								}
							if(!is.na(yvox) && !is.na(zvox){
								lambda[3,2]<-lambda[3,2]+(((vox+getPixels(img,x,y+1,z)[1])*(vox+getPixels(img,x,y,z+1)[1]))/(4*voxels*(subs-1))
								}
						}	
					}
				}
			}
		}
	lambda[2,1]<-lambda[1,2]
	lambda[3,1]<-lambda[1,3]
	lambda[2,3]<-lambda[3,2]

	fwhmx<-(4*log(2/lambda[1,1]))^(1/2)
	fwhmy<-(4*log(2/lambda[2,2]))^(1/2)
	fwhmz<-(4*log(2/lambda[3,3]))^(1/2)
	fwhm<-c(fwhmx,fwhmy,fwhmz)
##Tpdf<-probability density function of a t-distribution with v degress of freedom
	vintegral<-function(t){((((t^2)+subs-1)^2)/((v-1)*(v-2)))*((dt(t,df)^3)/(p(t)^2))}
	lambdav<-integrate(vintegral,Inf,Inf)
	lambdaz<-lambdav*lambdae
	return(fwhm,lambda)
	}
