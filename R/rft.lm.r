#' Fits a linear model according to random field theory
#' 
#' @param D image dimensions
#' @param x matrix or data frame of predictors
#' @param y observations (i.e. an image matrix)
#' @param conmat contrast matrix of contrasts * conditions (see example for more details)
#' @param mask object of type \code{antsImage}. Typically the same mask used to produce the image matrix.
#' @param statdir directory that statistical images are saved to (defaualt to "./")
#'
#'
#' @return a list of statistical images, coefficients, the FWHM estimate, and resel values.
#' @description
#'
#' This function is a wrapper for several functions necessary to performing a random field theory based
#' statistical analysis of images. A matrix ,or data frame, of predictors (\code{x}) is fitted to an 
#' image matrix (\code{y}) using a general linear model. The intercept is supressed and set to x=0. 
#' Therefore, \code{y} is scaled to zero. As recommended by ### et al. (1999) the residuals are used to 
#' estimate the full width at half-maxima (FWHM; or "smoothness") of the final statistical image. The 
#' the resolutions in voxel space (resels) are calculated using the estimated FWHM and the image search 
#' region (\code{mask}). The \code{T.imgs} output contains a list of t-statistical fields in image space
#' determined according to specified contrasts (\code{conmat}). 
#' 
#'
#' @References 
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worlsey K.J., et al. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Stefan J.K., (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' 
#' @Author Zachary P. Christensen
#' @kewords rft.pcluster, rft.ec, rft.resels
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
#' @examples
#' 
#' outimg1 <-makeImage(c(10,10,10), rnorm(1000))
#' maskimg <-getMask(outimg1)
#' fwhm <-est.Smooth(outimg1,maskimg)
#' # estimation of smoothness of overall sample images in a statistical model
#' outimg2 <-makeImage(c(10,10,10), rnorm(1000))
#' imat <-imageListToMatrix(list(outimg1, outimg2), maskimg)
#' var1 <-rnorm(2)
#' conmat <-matrix(c(1))
#' fit <-rft.fit(imat~var1-1, conmat, mask)
#' tmat <-fit$tfields
#' df <-fit$df
#' fwhm <-fit$fwhm
#' 
#'
#' thresh <-rft.thresh(3, timg, .05, 150, fwhm, mask, df, "T", "voxel")
#' results <-rft.results(3, thresh[1], 100, fwhm, timg, mask, df, "T", thresh, resels)
#'
#' @export rft.lm
rft.lm <-function(x, y, conmat, mask, tol=1e-07, statdir="./", findSmooth="TRUE", findResels="TRUE", concise="TRUE"){		
	z <-.lm.fit(x,y,tol=tol)
	p <-z$rank
	df <-c(n-p, p-1)
	b <-z$coefficients
	r <-z$residuals
	rss <-sum(r^2)
	mrss <-rss/df[1]
	p1 <-1L:p
	R <-chol2inv(z$qr[p1, p1, drop = FALSE])
	T.imgs <-list()
	for (i in 1:nrow(conmat)){  
		se <-t(as.matrix(sqrt(mrss *(conmat[i,] %*% R %*% conmat[i,]))))
		se[se==0] <-.01
		est <-conmat[i,] %*% b
		tstat <- est/se
		timg <-makeImage(mask,tstat)
		antsImageWrite(timg, file=paste(statdir, rownames(conmat)[i], ".nii.gz", sep=""))
		T.imgs <-lappend(T.imgs, timg)
		names(T.imgs)[i] <-rownames(conmat)[i]
		}
	if (findSmooth=="TRUE"){
		cat("Estimating FWHM. \n")
		fwhm <-estScaledSmooth(res,mask,df)
		}
	if (findResels=="TRUE"){
		cat("Calculating resels. \n")
		resels <-rft.resels(mask, fwhm)
		}
	if (concise=="TRUE"){
	ans <-list(design.matrix=x,
		T.imgs=T.imgs,
		conmat=conmat,
		df=df,
		fwhm=fwhm,
		resels=resels)
	}else {
	ans <-list(design.matrix=x,
		T.imgs=T.imgs,
		conmat=conmat,
		coefficients=coef,
		residuals=resid,
		df=df,
		fwhm=fwhm,
		resels=resels)
	}
	ans
	}
