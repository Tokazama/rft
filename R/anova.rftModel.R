#' Summarizing the F-statistic of rftModel fits
#' 
#' 
#' @param object object of class \code{rftGlm}
#' @param contrastMatrix optional matrix in which each row represents a unique contrast and each column the weights for each predictor
#' @param statdir optional directory where output is saved (if not specified images are not saved)
#' @param verbose enables verbose output
#' 
#' @return
#' \item{rftResults} {the outputs from \code{rftResults} for a particular contrast}
#' \item{contrastImage}
#' \item{dof} {degrees of freedom [effective degreesof freedom, residual degrees of freedom]}
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
#' @export anova.rftModel
anova.rftModel <- function(object, contrastMatrix, statdir = NULL, controls = list(), verbose = FALSE) {
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
  conNames <- rownames(contrastMatrix)
  if (is.null(conNames))
    conNames <- paste("contrast_", 1:ncon, sep = "")
  
  # compute contrasts----------------------------------------------------------
  RV <- object$modParams$R %*% object$modParams$V
  mrss <- object$rss / sum(diag(RV))
  conlist <- list()
  dof <- c(0, object$dof[2])
  for (i in 1:ncon) {
    C0 <- diag(object$dims$np) - c  %*% MASS::ginv(c)
    X0 <- X %*% C0
    R0 <- diag(object$dims$ns) - X0 %*% MASS::ginv(X0)
    M  <- R0 - R
    MV <- M %*% V
    traceMV <- sum(diag(MV))
    statvec <- ((t(object$coefficients) %*% t(X) %*% M %*% X %*% object$coefficients) / traceMV) / mrss
    # create image
    istat <- makeImage(object$dims$mask, statvec)
    if (!is.null(statdir)) {
      imgdir <- file.path(statdir, paste(conNames[i], "_", sep = ""))
      antsImageWrite(statimg, file = file.path(statdir, paste(statdir, conNames[i], ".nii.gz", sep = "")))
    } else
      imgdir <- NULL
    
    # calculate set/cluster/peak statistics------------------------------------
    dof[1] <- traceMV^2 / sum(diag(MV %*% MV))
    conlist[[i]] <- rftResults(istat, object$dims$resels, object$dims$fwhm, dof,
                               "F", rpvImage, controlvals$clusterThresh, controlvals$threshType,
                               controlvals$threshPval, controlvals$fdrThreshPval,
                               controlvals$conjuncImages, imgdir, verbose)
    
    conlist[[i]] <- structure(c(conlist[[i]], list(contrast = matrix(c, nrow = 1),
                                                   contrastImage = istat,
                                                   dof = dof)),
                              class = "anova.rftModel")
  }
  if (ncon > 1) {
    names(conlist) <- conNames
    structure(list(contrastResults = conlist,
                   fwhm = object$dims$fwhm,
                   resels = object$dims$resels,
                   voxels = object$dims$nv),
              class = "anova.rftModel.list")
  } else {
    c(conlist[[1]], list(fwhm = object$dims$fwhm,
                         resels = object$dims$resels,
                         voxels = object$dims$nv)
  }
}


#' @export print.anova.rftModel.list
print.anova.rftModel.list <- function(x) {
  ncon <- length(x$contrasts)
  
  for (i in 1:ncon) {
    cat(names(x$contrastResults)[i], "\n")
    cat("Contrast: ", x$contrastResults[i]$contrast, "\n")
    cat("Set-level: \n", sep = "")
    cat("  Clusters = ", x$contrastResults[[i]]$clusters, "\n")
    cat("  p-value = ", x$contrastResults[[i]]$set.level, "\n\n")
  }
  
  cat("Residual degrees of freedom = ", x$contrastResults[[1]]dof[2], "\n")
  cat("Cluster threshold = ", x$k, "\n")
  cat("FWHM = ", round(x$fwhm, 2), "\n")
  cat("Resels = ", round(x$resels), "\n")
  cat("Voxels = ", x$voxels, "\n")
}

#' @export print.anova.rftModel
print.anova.rftModel <- function(x) {
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
                   "Max F-stat"      = round(x$peak.level[, 4], 2),
                   "Z"               = round(x$peak.level[, 5], 2)))
  
  cat("Effective degrees of freedom = ", x$dof[1], "\n")
  cat("Residual degrees of freedom = ", x$dof[2], "\n")
  cat("Statistical threshold = ", round(x$threshold, 2), "\n")
  cat("Cluster threshold = ", x$contrastResults[[i]]$k, "\n")
  cat("FWHM = ", round(x$contrastResults[[i]]$fwhm, 2), "\n")
  cat("Resels = ", round(x$contrastResults[[i]]$resels), "\n")
  cat("Voxels = ", x$contrastResults[[i]]$voxels, "\n")
}
