
#' @export summary.rft
summary.rft <- function(object, contrast, conType, RPVImg = NULL, k = 1, threshType = "pRFT", pval = .05, pp = .001, n = 1, statdir, verbose = FALSE) {
  if (!is.null(RPVImg))
    rpvImage <- object$rpvImage
  if (missing(statdir))
    statdir <- object$statdir
  if (missing(contrast))
    contrast <- object$contrast
  if (missing(conType))
    conType = object$contrast.type
  if (is.null(contrast))
    contrast <- attr(object$X, "assign")
  if (class(contrast) == "numeric" | class(contrast) == "integer")
    contrast <- matrix(contrast, nrow = 1)
  # compute contrasts----------------------------------------------------------
  XX <- MASS::ginv(crossprod(object$dm))
  se <- sqrt((object$rss / object$df[2]) * as.numeric(contrast %*% tcrossprod(XX, contrast)))
  se[se == 0] <- .01 # control for NULLS that result from dividing by 0
  statmat <-(contrast %*% object$coefficients) / se
  if (conType =="F")
    statmat <- statmat^2
  istat <- makeImage(object$mask, statmat)
  
  # calculate set/cluster/peak statistics--------------------------------------
  ans <- rftResults(x = istat, resels = object$resels, fwhm = object$fwhm,
                    df = object$df, fieldType = conType, rpvImage = rpvImage,
                    k = k, threshType = threshType, pval = pval, pp = pp, n = n,
                    statdir = statdir, verbose = verbose) 
  ans <- list(set.level = ans$set.level,
              cluster.level = ans$cluster.level,
              peak.level = ans$peak.level,
              labeled.clusters = ans$labeled.clusters,
              threshold = ans$threshold,
              k = ans$k,
              clusterImage = ans$clusterImage,
              contrastImage = istat,
              contrast = contrast,
              contrast.type = conType,
              fwhm = object$fwhm,
              resels = object$resels,
              voxels = object$voxels,
              method = object$method)
  class(ans) <- "summary.rft"
  ans
  }

#' @export print.rft
print.rft <- function(x) {
  cat("Call: \n")
  cat(x$call, "\n\n")
  
  cat("Coefficients: \n")
  ncoef <- nrow(x$coefficients)
  for (i in 1:ncoef) {
    print(colnames(x$coef)[i], "\n")
  }
  
  cat("\nMethod = ", x$method, "\n", sep = "")
  cat("Voxels = ", x$voxels, "\n")
  cat("Resels = ", x$resels, "\n")
}




#' @export print.summary.rft
print.summary.rft <- function(x) {
  cat("Set-level: \n", sep = "")
  cat("Clusters = ", x$labeled.clusters, "\n")
  cat("p-value = ", x$set.level, "\n\n")
  
  cat("Cluster-Level: \n")
  print(data.frame("P-fwe"  = round(x$cluster.level[, 1], 3),
                   "P-fdr"  = round(x$cluster.level[, 2], 3),
                   "P"      = round(x$cluster.level[, 3], 3),
                   "Voxels" = x$cluster.level[, 4],
                   "x" = round(x$cluster.level[, 5]),
                   "y" = round(x$cluster.level[, 6]),
                   "z" = round(x$cluster.level[, 7])))

  cat("\nPeak-level: \n")
  cat("Field-type = ", x$contrast.type, "\n")
  print(data.frame("P-fwe" = round(x$peak.level[, 1], 3),
                   "P-fdr" = round(x$peak.level[, 2], 2),
                   "P"     = round(x$peak.level[, 3], 3),
                   "Max-Stat"      = round(x$peak.level[, 4], 2),
                   "Z"             = round(x$peak.level[, 5], 2)))
  
  cat("\nContrast = ", x$contrast, "\n")
  cat("Contrast type = ", x$contrast.type, "\n")
  cat("Statisticc threshold = ", round(x$threshold, 2), "\n")
  cat("Cluster threshold = ", x$k, "\n")
  cat("FWHM = ", round(x$fwhm, 2), "\n")
  cat("Resels = ", round(x$resels), "\n")
  cat("Voxels = ", x$voxels, "\n")
  cat("Method = ", x$method, "\n")
  }
