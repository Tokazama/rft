#' Estimates smoothness of an image
#' 
#' The partial derivatives of an image in x, y, and z directions are used
#' to create a covariance matrix which in turn is used to calculate the full-widths 
#' at half maxima (FWHM). The FWHM is equivalent to the estimated image smoothness.
#'
#' It is possible to use a single computed statistical parametric map image to
#' estimate the FWHM. However, it's recommended that FWHM estimates are obtained
#' from the residuals of statistical models (Stefan J.K et al., 1999). 
#'
#' @param img - object of class antsImage
#' @param mask - input mask, must match matrix
#' @return Outputs the estimated covariance matrix and fwhm in as an object of class list.
#' @references
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Stefan J.K., (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' @examples
#' # estimatation of a single images smoothness
#' outimg1 <-makeImage(c(10,10,10), rnorm(1000))
#' maskimg <-getMask(outimg1)
#' fwhm <-est.Smooth(outimg1,maskimg)
#' # estimation of smoothness of overall sample images in a statistical model
#' outimg2 <-makeImage(c(10,10,10), rnorm(1000))
#' imat <-imageListToMatrix(list(outimg1, outimg2), maskimg)
#' variable <-rnorm(2)
#' fit <-lm(imat ~ variable)
#' res <-residuals(fit)
#' DegFree <-fit$df.residual
#' meanres <-colMeans(res)
#' psd <-sqrt(colSums(res^2)/DegFree)
#' Zmat<-matrix(nrow=2,ncol=ncol(imat))
#' fwhm<-matrix(0L,nrow=1,ncol=3)
#' subs<-nrow(res)
#' for (i in 1:subs){
#' 	Zmat[i,]<-(res[i,]-Mres[1])/psd
#' 	img<-makeImage(mask,Zmat[i,])
#' 	smooth<-est.Smooth(img,mask)
#'	fwhm<-fwhm+smooth[[2]]
#' 	}
#' fwhm2<-sqrt(4*log(2)/(fwhm/(subs-1)))
#' 
#' @export est.Smooth
est.Smooth <-function(D, img, mask){
	if (D==2){
		dimx <-dim(img)[1]
		dimy <-dim(img)[2]
		d1 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		d2 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		m1 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		m2 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		maskar <-as.matrix(mask)
		voxels <-sum(maskar)
		imgar <-as.matrix(img)
		fwhm <-matrix(0L,nrow=1,ncol=3)
		lambda <-matrix(0L,nrow=1,ncol=3)

		#calculate partial derivatives y
		d1[2:(dimx+1), 2:(dimy+1)] <-imgar 
		d2[2:(dimx+1), 1:(dimy)] <-imgar
		m1[2:(dimx+1), 2:(dimy+1)] <-maskar
		m2[2:(dimx+1), 1:(dimy)] <-maskar
		m3 <- (m1 + m2) == 2
		Zyy <- ((d1 - d2)[m3 == 1])
		#variances of partial derivatives y
		lambda[1,1] <-sum(Zyy^2)/(voxels)

		#calculate partial derivatives y
		d1 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		d2 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		m1 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		m2 <-matrix(0, nrow=dimx+2, ncol=dimy+2)
		
		d1[2:(dimx+1), 2:(dimy+1)] <-imgar 
		d2[2:(dimx+1), 1:(dimy)] <-imgar
		m1[2:(dimx+1), 2:(dimy+1)] <-maskar
		m2[2:(dimx+1), 1:(dimy)] <-maskar
		m3 <- (m1 + m2) == 2
		Zyy <- ((d1 - d2)[m3 == 1])
		#variances of partial derivatives y
		lambda[1,2] <- sum(Zyy^2)/(voxels)

	}else{
		dimx <-dim(img)[1]
		dimy <-dim(img)[2]
		dimz <-dim(img)[3]
		d1 <-array(0, dim = dim(img) + 2)
		d2 <-array(0, dim = dim(img) + 2)
		m1 <-array(0, dim = dim(img) + 2)
		m2 <-array(0, dim = dim(img) + 2)
		maskar <-as.array(mask)
		voxels <-sum(maskar)
		imgar <-as.array(img)
		fwhm <-matrix(0L,nrow=1,ncol=3)
		lambda <-matrix(0L,nrow=1,ncol=3)

		#calculate partial derivatives x
		d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
		d2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
		m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
		m2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
		m3 <- (m1 + m2) == 2
		Zxx <- ((d1 - d2)[m3 == 1])
		#variances of partial derivatives x
		lambda[1,1] <- sum(Zxx^2)/(voxels)
	
		#calculate partial derivatives y
		d1 <-array(0, dim = dim(img) + 2)
		d2 <-array(0, dim = dim(img) + 2)
		m1 <-array(0, dim = dim(img) + 2)
		m2 <-array(0, dim = dim(img) + 2)

		d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar 
		d2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
		m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
		m2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
		m3 <- (m1 + m2) == 2
		Zyy <- ((d1 - d2)[m3 == 1])
		#variances of partial derivatives y
		lambda[1,2] <-sum(Zyy^2)/(voxels)

		#calculate partial derivatives z
		d1 <-array(0, dim = dim(img) + 2)
		d2 <-array(0, dim = dim(img) + 2)
		m1 <-array(0, dim = dim(img) + 2)
		m2 <-array(0, dim = dim(img) + 2)

		d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar 
		d2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
		m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
		m2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
		m3 <- (m1 + m2) == 2
		Zzz <- ((d1 - d2)[m3 == 1])
		#variances of partial derivatives z
		lambda[1,3] <-sum(Zzz^2)/(voxels)
		}

	#calculate fwhm from lambda
	fwhm <-sqrt(4*log(2)/lambda)
	smooth <-list("fwhm","lambda")
	smooth$fwhm <-fwhm
	smooth$lambda <-lambda
	return(smooth)
	}
