#
#' @name rft.lm
#' @title Fits design matrices
#' 
#' @param imat - Matrix of images of n x voxels
#' @param dm - Design matrix created by model.rft() function. With values dm["DesignMatrix", "DegreesOfFreedom", "D.M.Names"]
#' @param conmat - A contrast matrix of D.M.Names x contrast
#' @return Produces T-values, coefficients, and fwhm values.
#' @description
#'	Accepts a design matrix of 
#' \code{rft.lm} 
#'
#' @References Worlsey K.J., et al. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' @author Zachary P. Christensen
#' @kewords rft.pcluster, ants.ec
#' @examples
#' 
#' 
#'	
#' @export rft.lm
rft.lm <-function(imat,dm,conmat){
	if (missing(imat) | missing(dm) | missing(mask) | missing(conmat)){
		stop("Must specify imat dm and conmat")
	if (class(dm) !=matrix){
		stop("Design matrix must be of class matrix")
		}
	
	nsub <-nrow(imat)
	nvox <-ncol(imat)
	#Uses designm matrix to solve for T-values
	cols <-ncol(DM)
	U <-t(dm) %*% dm
	Usvd <- svd(U)
	if (is.complex(U)){
		Usvd$u <- Conj(Usvd$u)
		Positive <- Usvd$d > max(sqrt(.Machine$double.eps) * Usvd$d[1L], 0)
		}
	if (all(Positive)){ 
		UU <-Usvd$v %*% (1/Usvd$d * t(Usvd$u))
	}else if (!any(Positive)){
		UU <-array(0, dim(U)[2L:1L])
	}else{
		UU <-Usvd$v[, Positive, drop = FALSE] %*% ((1/Usvd$d[Positive]) * t(Usvd$u[, Positive, drop = FALSE]))
	}
	UY <-t(dm) %*% smat
	B <-UU %*% UY
	res <-smat - (dm %*% B)
	RSS <-colSums(res^2)
	MRSS <-RSS/DF
	tfields <-list()
	for (i in 1:ncol(conmat)){
		SE <-sqrt(MRSS *(conmat[,i] %*% UU %*% conmat[,i]))
		T<-(conmat[,i] %*% B)/SE
		tfields <-lappend(tfields,T,)
		}
	psd <-colMeans(sqrt(RSS))
	fwhm<-matrix(0L,nrow=1,ncol=3)
	Zmat <-matrix(nrow=nsub, ncol=nvox)
	cat("Estimating fwhm/smoothing")
	progress <- txtProgressBar(min = 0, max = nsub, style = 3)
	for (i in 1:nsub){
		Zmat<-(res[i,]-Mmat[1])/psd
		img<-makeImage(mask,Zmat[i,])
		smooth<-est.smooth(img,mask)
		fwhm<-fwhm+smooth[[2]]
		setTxtProgressBar(progress, i)
		}
	close(progress)
	fwhm2<-sqrt(4*log(2)/(fwhm/(nsub-1)))
	
	lmresults <-list(tfields,fwhm2)
	names(lmresults) <-c("tfields","fwhm")
	return(lmresults)
}
