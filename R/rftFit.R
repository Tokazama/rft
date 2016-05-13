#' Fit images using generalized linear models and random field theory
#' 
#' 
#' @param formula
#' @param mask
#' @param data
#' @param contrast
#' @param conType
#' @param statdir
#' @param sample
#' @param verbose
#' @param ...
#' 
#' @returnw
#' \item{call}
#' \item{coefficients} 
#' \item{residuals}
#' \item{rss}
#' \item{rpvImage} 
#' \item{fwhm}
#' \item{resels} 
#' \item{mask} 
#' \item{df}
#' \item{dm}
#' \item{contrast}
#' \item{conType}
#' \item{voxels}
#' \item{statdir}
#' 
#' @details
#' 
#' @author Zachary P. Christensen
#' 
#' 
#' 
#' 
#' 
#' 
#' @export summary.rftFit
rftFit <- function(formula, mask, data, contrast, conType, statdir = NULL, sample = NULL, verbose = TRUE, ...) {
  if (missing(data))
    data <- environment(formula)
  calll <- match.call()
  X <- model.matrix(formula, data = data)
  y <- model.response(model.frame(formula, data = data))
  
  v <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  edf <- p - 1
  rdf <- n - p
  
  b   <- matrix(0, nrow = p, ncol = v) # coefficients
  rownames(b) <- colnames(X)
  rss <- matrix(0, nrow = 1, ncol = v) # residual sum of squares
  r   <- matrix(0, nrow = n, ncol = v) # residuals
  
  # fit each voxel-------------------------------------------------------------
  if (verbose) {
    cat("Fitting model. \n")
    progress <- txtProgressBar(min = 0, max = v, style = 3)
  }
  for (i in 1:v) {
    fit <- glm.fit(X, y[, i], ...)
    b[, i]   <- fit$coefficients
    rss[, i] <- sum(fit$residuals^2)
    r[, i]   <- fit$residuals / sqrt(rss[, i] / rdf) # scale residuals with mean residual sum of squares
    if (verbose)
      setTxtProgressBar(progress, i)
  }
  if (verbose)
    close(progress)
  
  # estimate fwhm/resels-------------------------------------------------------
  if (!is.null(statdir)) {
    rl <- paste(statdir, "residual_", 1:n, ".nii.gz", sep = "")
    if (verbose)
      cat("Writing residual images to director. \n")
    for (i in 1:n) {
      antsImageWrite(makeImage(mask, r[i, ]), rl[i])
    }
    r <- rl
  }
  if (verbose)
    cat("Estimating FWHM. \n")
  smooth <- estSmooth(r, mask, rdf, scaleResid = FALSE, sample = sample, verbose = verbose)
  if (verbose)
    cat("Estimating resels. \n")
  v2r <- resels(mask, smooth$fwhm)
  
  if (missing(contrast))
    contrast <- NULL
  if (missing(conType))
    conType <- NULL
  # create rft object-----------------------------------------------------
  z <- structure(list(call = call,
                      coefficients = b,
                      residuals = r,
                      rss = rss,
                      RPVImage = smooth$RPVImage,
                      fwhm = smooth$fwhm,
                      resels = v2r,
                      mask = antsImageClone(mask),
                      df = c(edf, rdf),
                      dm = X,
                      contrast = contrast,
                      contrast.type = conType,
                      voxels = v,
                      statdir = statdir),
                 class = "rftFit")
  z
}



#' Summary of rftFit analysis
#' 
#' 
#' @param object object of class \code{rftFit}
#' @param contrast a matrix in which each row represents a unique contrast and each column the weights for each predictor
#' @param conType the type of contrast to compute ("F" or "T" (for "t"))
#' @param isotropic bool. should isotropy be assumed? (default TRUE)
#' @param k minimum desired cluster size (default = 1)
#' @param threshType a numeric value to threshold the statistical field or a character of the following methods:
#' \itemize{
#'	\item{cRFT: } {computes a threshold per expected cluster level probability }
#'	\item{pRFT: } {uses the mask and pval calculates the minimum statistical threshold}
#'	\item{cFDR: } {uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima}
#'	\item{pFDR: } {computes the fdr threshold for the entire field of voxels}
#' }
#' @param pval the p-value for estimating the threshold (default = .05)
#' @param pp the primary (initial) p-value for thresholding (only used for FDR methods; default = .001)
#' @param n number of images in conjunction
#' @param statdir directory where output is saved (if not specified images are not saved)
#' @param verbose enables verbose output
#' 
#' @return
#' \item{call} {original call to \code{rftFit}}
#' \item{contrast.matrix} {matrix used for contrasts}
#' \item{conType} {the type of statistical field for each contrast}
#' \item{contrastResults} {Contains the resulting image for each contrast along with all the results produce by \code{rftResults} for that particular contrast.}
#' \item{fwhm} {full width at half maxima}
#' \item{resels} {resolution elements}
#' \item{voxels} {number of active voxels (number of voxels in mask)}
#' 
#' @details
#' 
#' @author Zachary P. Christensen
#' 
#' 
#' 
#' 
#' 
#' 
#' @export summary.rftFit
summary.rft <- function(object, contrast, conType, isotropic = TRUE, k = 1,
                        threshType = "pRFT", pval = .05, pp = .001, n = 1,
                        statdir, verbose = FALSE) {
  if (isotropic == FALSE)
    rpvImage <- object$rpvImage
  if (missing(statdir))
    statdir <- object$statdir
  if (missing(contrast))
    contrast <- object$contrast
  else {
    if (class(contrast) != "matrix" | class(contrast) != "dsCMatrix" | class(contrast) != "data.frame")
      contrast <- as.matrix(contrast, nrow = 1)
    ncon <- nrow(contrast)
    if (ncol(object$dm) != ncol(contrast))
      stop("Each contrast length must be equal to the number of columns in the design matrix, including any intercepts. \n")
  }
  if (missing(conType))
    conType = object$conType
  if (is.null(contrast))
    contrast <- attr(object$X, "assign")
  if (class(contrast) == "numeric" | class(contrast) == "integer")
    contrast <- matrix(contrast, nrow = 1)
  if (is.null(rownames(contrast))) {
    conNames <- paste("contrastImage", 1:ncon, sep = "")
    rownames(contrast) <- conNames
  }
  else
    conNames <- rownames(contrast)
  if (ncon > 1) {
    if (length(k) != ncon)
      k <- rep(k[1], ncon)
    if (length(pval) != ncon)
      pval <- rep(pval[1], ncon)
    if (length(pp) != ncon)
      pp <- rep(pp[1], ncon)
  }
  # compute contrasts----------------------------------------------------------
  XX <- MASS::ginv(crossprod(object$dm))
  cons <- list()
  for (i in 1:ncon) {
    c <- matrix(contrast[i, ], nrow = 1)
    se <- sqrt((object$rss / object$df[2]) * as.numeric(c %*% tcrossprod(XX, c)))
    se[se == 0] <- .01 # control for NULLS that result from dividing by 0
    statmat <-(c %*% object$coefficients) / se
    if (conType =="F")
      statmat <- statmat^2
    istat <- makeImage(object$mask, statmat)
    if (!is.null(statdir))
      antsImageWrite(statimg, file = paste(statdir, conNames[i], ".nii.gz", sep=""))
  
    # calculate set/cluster/peak statistics--------------------------------------
    ans <- rftResults(ilist[[i]], resels = object$resels, fwhm = object$fwhm,
                      df = object$df, fieldType = conType[i], rpvImage = rpvImage,
                      k = k[i], threshType = threshType[i], pval = pval[i], pp = pp[i],
                      n = n, statdir = statdir, verbose = verbose)
    
    cons <- lappend(cons, list(set.level = ans$set.level,
                          cluster.level = ans$cluster.level,
                          peak.level = ans$peak.level,
                          clusters = ans$clusters,
                          threshold = ans$threshold,
                          k = ans$k,
                          clusterImage = ans$clusterImage,
                          contrastImage = istat))
  }
  
  names(ilist) <- conNames
  colnames(contrast) <- rownames(object$coefficients)
  structur(list(call = object$call,
                contrast.matrix = contrast,
                conType = conType,
                contrastResults = cons,
                fwhm = object$fwhm,
                resels = object$resels,
                voxels = object$voxels),
           class = "summary.rftFit")
}


#' @export print.rftFit
print.rft <- function(x) {
  cat("Call: \n")
  cat(x$call, "\n\n")
  
  cat("Coefficients: \n")
  ncoef <- nrow(x$coefficients)
  for (i in 1:ncoef) {
    print(colnames(x$coefficients)[i], "\n")
  }
  
  cat("Voxels = ", x$voxels, "\n")
  cat("Resels = ", x$resels, "\n")
}

#' @export print.summary.rftFit
print.summary.rftFit <- function(x, display) {
  ncon <- length(x$contrasts)
  if (missing(display)) {
    if (ncon > 1)
      display <- 2
    else
      display <- 1
  }
  cat("Contrast Matrix: \n")
  conmat <- cbind(object$conType, object$contrast.matrix)
  colnames(conmat)[1] <- "Contrast-type"
  print(as.data.frame(conmat), "\n")
  
  for (i in 1:ncon) {
    cat(names(x$contrastResults[[i]]), "\n")
    cat("Set-level: \n", sep = "")
    cat("Clusters = ", x$contrastResults[[i]]$clusters, "\n")
    cat("p-value = ", x$contrastResults[[i]]$set.level, "\n\n")
    if (display == 1) {
      cat("Cluster-Level: \n")
      print(data.frame("P-fwe"  = round(x$contrastResults[[i]]$cluster.level[, 1], 3),
                       "P-fdr"  = round(x$contrastResults[[i]]$cluster.level[, 2], 3),
                       "P"      = round(x$contrastResults[[i]]$cluster.level[, 3], 3),
                       "Voxels" = x$contrastResults[[i]]$cluster.level[, 4],
                       "x" = round(x$contrastResults[[i]]$cluster.level[, 5]),
                       "y" = round(x$contrastResults[[i]]$cluster.level[, 6]),
                       "z" = round(x$contrastResults[[i]]$cluster.level[, 7])))
      
      cat("\nPeak-level: \n")
      cat("Field-type = ", x$contrastResults[[i]]$conType, "\n")
      print(data.frame("P-fwe" = round(x$contrastResults[[i]]$peak.level[, 1], 3),
                       "P-fdr" = round(x$contrastResults[[i]]$peak.level[, 2], 2),
                       "P"     = round(x$contrastResults[[i]]$peak.level[, 3], 3),
                       "Max-Stat"      = round(x$contrastResults[[i]]$peak.level[, 4], 2),
                       "Z"             = round(x$contrastResults[[i]]$peak.level[, 5], 2)))
      
      cat("Statistical threshold = ", round(x$contrastResults[[i]]$threshold, 2), "\n")
      cat("Cluster threshold = ", x$contrastResults[[i]]$k, "\n")
    }
    cat("FWHM = ", round(x$contrastResults[[i]]$fwhm, 2), "\n")
    cat("Resels = ", round(x$contrastResults[[i]]$resels), "\n")
    cat("Voxels = ", x$contrastResults[[i]]$voxels, "\n")
  }
}
