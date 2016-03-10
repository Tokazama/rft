#' Estimates smoothness of an image
#' 
#' The partial derivatives of an image in x, y, and z directions are used
#' to create a covariance matrix which in turn is used to calculate the full-widths
#' at half maxima (FWHM). The FWHM is equivalent to the estimated image smoothness.
#' 
#' The resels per voxel image (\code{RPVImg}) represents the estimated resel at each 
#' individual voxel throughout the search region. This may be used in place of volumetric
#' measurements (or sum voxel measurements) when estimating the p-value of a cluster
#' using \code{rftPval}. The intent behind using the RPV image to estimate cluster level 
#' statistics is to offset the natural probability of obtaining significant clusters solely by
#" chance in very smooth regions at low thresholds.
#'
#' It is possible to use a single statistical field image to estimate the FWHM. However, it's
#' recommended that FWHM estimates are obtained from the scaled residuals of statistical models
#' (Stefan J.K et al., 1999). Therefore, this function is optimized to estimate the pooled 
#' smoothness of the residual images from a fitted model. 
#' 
#' A scaling factor is used to correct for differences when using the \code{sample} option.
#' Scaling isn't effective when the number of images is very low and typically results in an 
#' overestimation of the the FWHM. If only one image or numeric vector is entered then the 
#' scaling factor is not used. If a numeric vector is entered the \code{imageMake} function is 
#' used to prepare it for smoothness estimation (see Worsley et al., 1999).
#'
#' @param x object of class antsImage
#' @param mask input mask, must match matrix
#' @param rdf residual degrees of freedom
#' @param scaleResid scale residuals?
#' @param makeRPV make resels per voxel (RPV) image?
#' @param sample number of images to use for estimating smoothing (default uses all images)
#' @param verbose enables verbose output
#' 
#' @return Outputs the estimated FWHM and RPV image (if \code{makeRPV = TRUE})
#' 
#' @references
#' Hayasaka (2004) Nonstationary cluster-size inference with random field and permutation methods.
#' Worsley K.J. (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worsley K.J. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Worsley K.J. (1999) Detecting Changes in Nonisotropic Images
#' Stefan J.K. (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' 
#' @Author Zachary P. Christensen
#' 
#' @seealso resels
#' 
#' @note: function currently in beta phase
#' @examples
#'
#' # estimatation of a single images smoothness
#' outimg1 <- makeImage(c(10, 10, 10), rnorm(1000))
#' maskimg <- getMask(outimg1)
#' myfwhm1 <- estSmooth(outimg1, maskimg)
#' 
#' # estimation of smoothness of overall sample images in a statistical model
#' outimg2 <- makeImage(c(10,10,10), rnorm(1000))
#' imat <- imageListToMatrix(list(outimg1, outimg2), maskimg)
#' variable <- rnorm(2)
#' fit <- lm(imat ~ variable)
#' myfwhm2 <- estSmooth(residuals(fit), maskimg)
#' 
#' @export estSmooth
estSmooth <- function(x, mask, rdf, scaleResid = TRUE, makeRPV = FALSE, sample, verbose = FALSE) {
    D <- mask@dimension

    dimx  <- 1:dim(mask)[1]
    dimx1 <- 2:(dim(mask)[1] + 1)
    dimx2 <- 3:(dim(mask)[1] + 2)
    if (D > 1) {
        dimy  <- 1:dim(mask)[2]
        dimy1 <- 2:(dim(mask)[2] + 1)
        dimy2 <- 3:(dim(mask)[2] + 2)
    }
    if (D > 2) {
        dimz  <- 1:dim(mask)[3]
        dimz1 <- 2:(dim(mask)[3] + 1)
        dimz2 <- 3:(dim(mask)[3] + 2)
    }

    if (class(x) == "antsImage") {
        if (D > 1)
            imgar <- as.matrix(x)
        if (D > 2)
            imgar <- as.array(x)
        scale <- 1
        n <- 1
        mrss <- 1
    } else if (class(x) == "numeric") {
        x <- matrix(x, nrow = 1)
        scale <- 1
        n <- 1
        mrss <- 1
    } else if (class(x) == "matrix") {
        mrss <- sqrt(colSums(x ^ 2) / rdf) # mean residual sum of squares for standardizing images later on
        if (missing(sample)) {
            nfull <- nrow(x) # original number of images (rows)
        } else {
            nfull <- nrow(x)
            rsamples <- sample(nrow(x), sample)
            x <- x[rsamples,]
        }
        n <- nrow(x) # number of images in sample (rows)
        scale <- (nfull / df) * (1 / n)
    }

    if (makeRPV == TRUE)
        if (D > 1)
            Vxy <- matrix(0, nrow = dim(mask)[1], ncol = dim(mask)[2])
    else if (D > 2)
        Vxy <- Vxz <- Vyz <- array(0, dim = dim(mask))
    if (D > 1) {
        maskar <- as.array(mask)
        Vxx <- Vyy <- matrix(0, nrow = dim(mask)[1], ncol = dim(mask)[2])
        m <- d <- matrix(0, nrow = dim(mask)[1] + 2, ncol = dim(mask)[2] + 2)
        xyzm <- m[dimx, dimy1] + m[dimx1, dimy1] + m[dimx2, dimy1] +
                m[dimx1, dimy] + m[dimx1, dimy1] + m[dimx1, dimy2]
        xyzm[xyzm != 6] <- 0
        xyzm[xyzm == 6] <- 1
        imgar <- maskar
    } else if (D > 2) {
        maskar <- as.array(mask)
        Vxx <- Vyy <- Vzz <- array(0, dim = dim(mask))
        m <- d <- array(0, dim = dim(mask) + 2)
        m[dimx1, dimy1, dimz1] <- maskar
        xyzm <- m[dimx, dimy1, dimz1] + m[dimx1, dimy1, dimz1] + m[dimx2, dimy1, dimz1] +
                m[dimx1, dimy, dimz1] + m[dimx1, dimy1, dimz1] + m[dimx1, dimy2, dimz1] +
                m[dimx1, dimy1, dimz] + m[dimx1, dimy1, dimz1] + m[dimx1, dimy1, dimz2]
        xyzm[xyzm != 9] <- 0
        xyzm[xyzm == 9] <- 1
        imgar <- maskar
    }

    nvox <- sum(xyzm)

    if (verbose)
        progress <- txtProgressBar(min = 0, max = n, style = 3)
    for (i in 1:n) {
        if (class(x) == "matrix")
            imgar[maskar == 1] <- x[i,] / mrss
        if (D == 1) {
            d[dimx1] <- imgar
            dx <- (imgar - d[dimx])
        } else if (D == 2) {
            d[dimx1, dimy1] <- imgar
            dx <- (imgar - d[dimx, dimy1])
            dy <- (imgar - d[dimx1, dimy])
        } else if (D > 3) {
            d[dimx1, dimy1, dimz1] <- imgar
            dx <- (imgar - d[dimx, dimy1, dimz1])
            dy <- (imgar - d[dimx1, dimy, dimz1])
            dz <- (imgar - d[dimx1, dimy1, dimz])
        }

        Vxx <- Vxx + (dx * dx)

        if (D > 1)
            Vyy <- Vyy + (dy * dy)
        if (D > 2)
            Vzz <- Vzz + (dz * dz)
        if (makeRPV == "TRUE") {
            if (D > 1)
                Vxy <- Vxy + (dx * dy)
            if (D > 2) {
                Vxz <- Vxz + (dx * dz)
                Vyz <- Vyz + (dy * dz)
            }
        }
        if (verbose)
            setTxtProgressBar(progress, i)
    }
    if (verbose)
        close(progress)

    Vxx <- Vxx * scale
    if (D > 1)
        Vyy <- Vyy * scale
    if (D > 2)
        Vzz <- Vzz * scale
    if (makeRPV == "TRUE") {
        if (D > 1)
            Vxy <- Vxy * scale
        if (D > 2) {
            Vxz <- Vxz * scale
            Vyz <- Vyz * scale
        }
        rpv <- Vxx * Vyy * Vzz + Vxy * Vyz * Vxz * 2 - Vyz * Vyz * Vxx - Vxy * Vxy * Vzz - Vxz * Vxz * Vyy
        rpv[rpv < 0] <- 0
        rpv <- sqrt(rpv / (4 * log(2)) ^ D)
        rpv <- rpv * maskar
        RPVImg <- makeImage(mask, rpv)
    }
    if (D == 2)
        xyz <- cbind(matrix((Vxx * xyzm), ncol = 1), matrix((Vyy * xyzm), ncol = 1))
    else if (D > 2)
        xyz <- cbind((Vxx * xyzm), (Vyy * xyzm), (Vzz * xyzm))
    xyz <- sqrt(xyz / (4 * log(2)))
    xyz <- colSums(xyz) / nvox
    if (makeRPV == "TRUE") {
        rpv <- sum(rpv) / nvox
        R <- rpv ^ (1 / D) * (xyz / prod(xyz) ^ (1 / D))
        fwhm <- 1 / R
        list(fwhm = fwhm, RPVImg = RPVImg)
    } else if (classval == 1 | classval == 2) {
        fwhm <- 1 / xyz
        list(fwhm = fwhm, RPVImg = NULL)
    }
}