#
#' rft.lm
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
rft.lm <-function(formula, conmat, mask, test="FALSE"){	
	dm <-model.matrix(formula)
	nsub <-nrow(dm)
	nvar <-ncol(dm)
	if (test=="TRUE"){
		z <-list(design.matrix,contrast.matrix)
		z$design.matrix <-dm
		conmat <-matrix(0L,nrow=nvar, ncol=nvar)
		rownames(conmat) <-colnames(dm)
		z$contrast.matrix <-conmat
	}else{
		dm <-cbind(dm,1)
		conmat <-rbind(conmat,0)
		DF <-nsub-nvar
		imat <-model.response(model.frame(formula))
		nvox <-ncol(imat)
		#Uses designm matrix to solve for T-values
		U <-t(dm) %*% dm
		Usvd <- svd(U)
		Usvd$u <- Conj(Usvd$u)
		Positive <- Usvd$d > max(sqrt(.Machine$double.eps) * Usvd$d[1L], 0)
		if (all(Positive)){ 
			UU <-Usvd$v %*% (1/Usvd$d * t(Usvd$u))
		}else if (!any(Positive)){
			UU <-array(0, dim(U)[2L:1L])
		}else{
			UU <-Usvd$v[, Positive, drop = FALSE] %*% ((1/Usvd$d[Positive]) * t(Usvd$u[, Positive, drop = FALSE]))
		}
		UY <-t(dm) %*% imat
		B <-UU %*% UY
		cat("Calculating residuals
		",sep="")
		res <-imat - (dm %*% B)
		RSS <-colSums(res^2)
		MRSS <-RSS/DF
		tfields <-matrix(nrow=ncol(conmat),ncol=nvox)
		for (i in 1:ncol(conmat)){
			SE <-t(as.matrix(sqrt(MRSS *(conmat[,i] %*% UU %*% conmat[,i]))))
			SE[SE==0] <-.01
			tfields[i,] <-(conmat[,i] %*% B)/SE
			}
		rownames(tfields) <-colnames(conmat)
		psd <-colMeans(sqrt(as.matrix(RSS)))
		fwhm<-matrix(0L,nrow=1,ncol=3)
		Mmat <-colMeans(res)
		Zmat <-matrix(nrow=nsub, ncol=nvox)
		cat("Estimating fwhm/smoothing
		",sep="")
		progress <- txtProgressBar(min = 0, max = nsub, style = 3)
		for (i in 1:nsub){
			Zmat[i,]<-(res[i,]-Mmat[1])/psd
			img<-makeImage(mask,Zmat[i,])
			smooth<-est.Smooth(img,mask)
			fwhm<-fwhm+smooth[[2]]
			setTxtProgressBar(progress, i)
			}
		close(progress)
		fwhm2<-sqrt(4*log(2)/(fwhm/DF))
		cat("Calculating resels
		",sep="")
		resels <-rft.resel(mask, fwhm2)
		z <-list("design.matrix", "tfields", "coefficients", "df", "fwhm", "residuals", "contrast.matrix", "resels")
		z$design.matrix <-dm
		z$tfields <-tfields
		z$coefficients <-B
		z$df <-DF
		z$fwhm <-fwhm2
		z$residuals <-res
		z$contrast.matrix <-conmat
		z$resels <-resels
		}
	z
	return(z)
}
