#' Estimate image smoothness
#'
#' Estimates smoothness of a single image or image matrix
#'
#' @param x May be an image of class "antsImage" or an image matrix.
#' @param mask Input mask, must match matrix.
#' @param rdf Residual degrees of freedom.
#' @param scaleResid If \code{TRUE} residuals are scaled.
#' @param sample Number of images to use for estimating smoothing (default uses all images).
#' 
#' @return Outputs the estimated FWHM and resel per voxel image.
#' 
#' @details
#' The partial derivatives of an image in x, y, and z directions are used to
#' create a covariance matrix which in turn is used to calculate the 
#' full-widths at half maxima (FWHM). The FWHM is equivalent to the estimated
#' image smoothness.
#' 
#' The resels per voxel image (\code{rpvImage}) represents the estimated resel at
#' each individual voxel throughout the search region. This may be used in 
#' place of volumetric measurements (or sum voxel measurements) when estimating
#' the p-value of a cluster using \code{\link{rftPval}}. The intent behind using the 
#' RPV image to estimate cluster level statistics is to offset the natural
#' probability of obtaining significant clusters solely by chance in very 
#' smooth regions at low thresholds.
#'
#' It is possible to use a single statistical field image to estimate the FWHM. 
#' However, it's recommended that FWHM estimates are obtained from the scaled 
#' residuals of statistical models (Stefan J.K et al., 1999). Therefore, this 
#' function is optimized to estimate the pooled smoothness of the residual 
#' images from a fitted model.
#' 
#' A scaling factor is used to correct for differences when using the 
#' \code{sample} option. Scaling isn't effective when the number of images is 
#' very low and typically results in an overestimation of the the FWHM. If only
#' one image or numeric vector is entered then the scaling factor is not used. 
#' If a numeric vector is entered the \code{imageMake} function is used to 
#' prepare it for smoothness estimation (see Worsley et al., 1999).
#' 
#' Any NA values in \code{object} will be set to zero.
#' 
#' @references
#' Hayasaka (2004) Nonstationary cluster-size inference with random field and permutation methods.
#' 
#' Worsley K.J. (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' 
#' Worsley K.J. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' 
#' Worsley K.J. (1999) Detecting Changes in Nonisotropic Images
#' 
#' Stefan J.K. (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' 
#' @author Zachary P. Christensen
#' 
#' @seealso \code{\link{resels}}
#' 
#' @examples
#' # estimate individual image
#' mnit1 <- antsImageRead(getANTsRData('mni'))
#' mask <- getMask(mnit1)
#' fwhm1 <- estSmooth(mnit1, mask)
#' 
#' @export estSmooth
estSmooth <- function(x, mask, rdf, scaleResid = TRUE, sample = NULL) {
  DIM <- dim(mask)
  nfull <- nrow(x)
  if (!is.null(sample))
   s <- sample(1:nfull, sample)
  else
    s <- nfull
  smooth <- .Call("estSmooth", x[s, ], as.array(mask), rdf, nfull, DIM, scaleResid, package = "ANTsR")
  smooth$rpvImage <- makeImage(mask, smooth$rpvImage)
  smooth
  }
