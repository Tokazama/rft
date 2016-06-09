#' Summarizing the t-statistic of rftModel fits
#' 
#' 
#' @param object object of class \code{rftGlm}
#' @param contrastMatrix optional matrix in which each row represents a unique contrast and each column the weights for each predictor
#' @param statdir optional directory where output is saved (if not specified images are not saved)
#' @param verbose enables verbose output
#' 
#' @return
#' \item{call} {original call to \code{rftFit}}
#' \item{contrastmatrix} {matrix used for contrasts}
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
#' @export summary.rftGlm
summary.rftModel <- function(object, contrastMatrix, statdir = NULL, controls = list(), verbose = FALSE) {
  # check arguments------------------------------------------------------------
  controlvals <- rftControl()
  if (!missing(control))
    controlvals[names(control)] <- control
  if (controlvals$isotropic)
    rpvImage <- object$dims$rpvImage
  if (missing(contrastMatrix))
    contrastMatrix <- object$contrastMatrix
  else {
    condim <- dim(contrastMatrix)
    if (is.null(condim))
      contrastMatrix <- matrix(contrastMatrix, 1)
  }
  ncon <- nrow(contrastMatrix)
  
  if (ncol(contrastMatrix) != object$dims$np)
    stop("Each contrast length must be equal to the number of columns in the design matrix, including any intercepts. \n")
  if (is.null(rownames(contrastMatrix))) {
    conNames <- paste("contrastImage", 1:ncon, sep = "")
    rownames(contrastMatrix) <- conNames
  }
  else
    conNames <- rownames(contrastMatrix)
  if (is.null(rownames(contrastMatrix)))
    conNames <- paste("contrastImage", 1:ncon, sep = "")
  else
    conNames <- rownames(contrast)
  
  # compute contrasts----------------------------------------------------------
  XX <- MASS::ginv(crossprod(object$X))
  RV <- object$R %*% object$u$V
  mrss <- object$rss / sum(diag(RV))
  conlist <- list()
  dof <- c(0, object$dof[2])
  for (i in 1:ncon) {
    c <- matrix(contrast[i, ], ncol = 1)
    se <- sqrt(mrss * (crossprod(c, XX) %*% c))
    statvec <- crossprod(c, object$beta) / se
    
    # create image
    istat <- makeImage(object$dims$mask, statvec)
    if (!is.null(statdir))
      antsImageWrite(statimg, file = file.path(statdir, paste(conNames[i], ".nii.gz", sep = "")))
    
    # calculate set/cluster/peak statistics--------------------------------------
    conlist[[i]] <- rftResults(istat, object$dims$resels, object$dims$fwhm, dof,
                               "T", rpvImage, controlvals$clusterThresh, controlvals$threshType,
                               controlvals$threshPval, controlvals$fdrThreshPval,
                               controlvals$conjuncImages, statdir, verbose)
    conlist[[i]] <- structure(c(conlist[[i]], list(conNames[i] = istat, dof = dof)), class = "summary.rftModel")
  }
  if (ncon > 1) {
    names(conlist) <- conNames
    structure(list(contrastResults = conlist,
                   fwhm = object$dims$fwhm,
                   resels = object$dims$resels,
                   voxels = object$dims$nv),
              class = "summary.rftModel.list")
  } else {
    c(conlist[[1]], list(contrastMatrix = contrastMatrix,
                         fwhm = object$dims$fwhm,
                         resels = object$dims$resels,
                         voxels = object$dims$nv)
  }
}

#' @export print.summary.rftModel.list
print.summary.rftModel.list <- function(x) {
  ncon <- length(x$contrasts)
  cat("Contrast Matrix: \n")
  print(object$contrastMatrix, "\n")
  
  for (i in 1:ncon) {
    cat(names(x$contrastResults)[i], "\n")
    cat("Contrast: ", x$contrastResults[i]$contrast, "\nj")
    cat("Set-level: \n", sep = "")
    cat("  Clusters = ", x$contrastResults[[i]]$clusters, "\n")
    cat("  p-value = ", x$contrastResults[[i]]$set.level, "\n\n")
    }
  
  cat("Degrees of freedom = ", x$constrastResults[[1]]dof[2], "\n")
  cat("Cluster threshold = ", x$k, "\n")
  cat("FWHM = ", round(x$fwhm, 2), "\n")
  cat("Resels = ", round(x$resels), "\n")
  cat("Voxels = ", x$voxels, "\n")
}

#' @export print.summary.rftModel
print.summary.rftModel <- function(x) {
  cat("Set-level: \n", sep = "")
  cat("Clusters = ", x$clusters, "\n")
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
  print(data.frame("P-fwe" = round(x$peak.level[, 1], 3),
                   "P-fdr" = round(x$peak.level[, 2], 2),
                   "P"     = round(x$peak.level[, 3], 3),
                   "Max t-stat"      = round(x$peak.level[, 4], 2),
                   "Z"               = round(x$peak.level[, 5], 2)))
  
  cat("Residual degrees of freedom = ", x$dof[2], "\n")
  cat("Statistical threshold = ", round(x$threshold, 2), "\n")
  cat("Cluster threshold = ", x$k, "\n")
  cat("FWHM = ", round(x$fwhm, 2), "\n")
  cat("Resels = ", round(x$resels), "\n")
  cat("Voxels = ", x$voxels, "\n")
}
