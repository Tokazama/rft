#' Fits a linear model according to random field theory
#' 
#' @param x matrix or data frame of predictors
#' @param y observations (i.e. an image matrix)
#' @param conmat contrast matrix of contrasts by conditions (see example for more details)
#' @param mask object of type \code{antsImage}. Typically the same mask used to produce the image matrix.
#' @param tol
#' @param statdir directory that statistical images are saved to (defaualt to "./")
#' @param findSmooth if \cod{findSmooth="FALSE"} then FWHM and resels aren't calculated (default is \code{"TRUE"}
#' @param concise if \code{concise="FALSE"} then residuals and coefficients are returned (default is \code{"TRUE"})
#' @param resSample how many of the residual images will be used to estimate FWHM and resels (default is all images)
#'
#' @return A list comprising:
#' \describe{
#' \item{conmat}{contrast matrix used to create \code{Timgs}}
#' \item{Designmatrix}{the design matrix used to fit the model}
#' \item{df}{degrees of freedom}
#' \item{residuals}{image matrix of residuals}
#' \item{coefficients}{image matrix of coefficients}
#' \item{Timgs}{list of the t-field images } 
#' \item{fwhm}{the estimated FWHM}
#' \item{RPVImg}{resels per voxel image}
#' \item{resels}{estimated resels for the search space}
#' }
#' @description
#'
#'
#' This function is a wrapper for several functions necessary to performing a random field theory based
#' statistical analysis of images. A matrix ,or data frame, of predictors (\code{x}) is fitted to an 
#' image matrix (\code{y}) using \code{.lm.fit}. The images are each mean scaled to zero prior to
#' fitting (Friston K.J., 1995). As recommended by Stefan et al. (1999), the residuals are used to 
#' estimate the full width at half-maxima (FWHM; or "smoothness") of the final statistical image. The
#' residuals are standardized using the mean sum of squares of the residuals.
#' 
#' Standardized Residuals = residuals/sqrt(MRSS/Error of Degrees of Freedom)
#' 
#' The the resolutions in voxel space (resels) are calculated using the estimated FWHM and the image search 
#' region (\code{mask}). The \code{T.imgs} output contains a list of t-statistical fields in image space
#' determined according to specified contrasts (\code{conmat}).
#' 
#'
#'
#' @References 
#' Friston K.J., (1995) Statistical Parametric Maps in Functional Imaging: A General Linear Approach 
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
rft.lm <-function(x, y, conmat, mask, tol=1e-07, statdir="./", findSmooth="TRUE", concise="TRUE",resSample){		
  z <-.lm.fit(x,y,tol=tol)
  p <-z$rank
  df <-c(p-1,nrow(z$residuals)-p)
  b <-z$coefficients
  r <-z$residuals
  mrss <-colSums(r^2)/df[2]
  p1 <-1L:p
  R <-chol2inv(z$qr[p1, p1, drop = FALSE])  
  T.imgs <-list()
  if (concise=="TRUE"){
    ans <-list(conmat=conmat,DesignMatrix=x,df=df)
  }else {
    ans <-list(conmat=conmat,DesignMatrix=x,df=df,residuals=r,coefficients=b)
  }
  if (class(conmat)=="numeric"){
    conmat <-matrix(conmat, nrow=1)
  }
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
  ans <-lappend(ans,T.imgs)
  names(ans)[length(ans)] <-"Timgs"
  if (findSmooth=="TRUE"){
    cat("Estimating FWHM. \n")
    sr <-r/sqrt(mrss) #standardize residuals
    if (missing(resSample)){
      fwhm <-estSmooth(sr,mask,df)
    }else{
      fwhm <-estSmooth(sr,mask,df,resSample)
    }
    ans <-lappend(ans,fwhm$fwhm)
    names(ans)[length(ans)] <-"fwhm"
    ans <-lappend(ans,fwhm$RPVImg)
    names(ans)[length(ans)] <-"RPVImg"
    
    cat("Calculating resels. \n")
    resels <-rft.resels(mask, fwhm$fwhm)
    ans <-lappend(ans,resels)
    names(ans)[length(ans)] <-"resels"
  }
  ans
  }
