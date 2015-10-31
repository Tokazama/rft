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
	partial.derivative<-function(img,x,y,z,lambda){
		vox<-getPixels(img,x,y,z)[1]
		if(!is.na(vox)){
			X<-getpixels(img,x+1,y,z)[1]
			}
		if(!is.na(xvox)){
			lambda[1,1]<-((X-vox)^2)/(N*(n-1))
			}
		Y<-getpixels(img,x,y+1,z)[1]
		if(!is.na(Y)){
			lambda[2,2]<-((Y-vox)^2)/(N*(n-1))
			}
		Z<-getpixels(img,x,y,z+1)[1]
		if(!is.na(Z)){
			lambda[3,3]<-((Z-vox)^2)/(N*(n-1))
			}
		XY<-getPixels(img,x+1,y+1,z)[1]
		if(!is.na(X) && !is.na(Y) && !is.na(XY)){
			lambda[1,2]<-(((X-vox)+(Y-XY))*((Y-vox)+(X-XY)))/(4*N*(n-1))
			}
		XZ<-getPixels(img,x+1,y,z+1)[1]
		if(!is.na(X) && !is.na(Z) && !is.na(XZ)){
			lambda[1,3]<-(((X-vox)+(Z-XZ))*((Z-vox)+(X-XZ)))/(4*N*(n-1))
			}
		YZ<-getPixels(img,x,y+1,z+1)[1]
		if(!is.na(Y) && !is.na(Z) && !is.na(YZ)){
			lambda[2,3]<-(((Y-vox)+(Y-YZ))*((Z-vox)+(Y-YZ)))/(4*N*(n-1))
			}
		}

	for (i in 1:subs){
		img<-makeImage(mask,Zmat[i,]
		for (x in 1:(dimx)){
			for (y in 1:(dimy)){
				for (z in 1:(dimz)){
					lambda<-partial.derivative(img,x,y,z,lambda)
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
