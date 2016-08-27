# to do:
# multiple modality eigenanatomy stuff
# impute for ilm

#' iData Model Formulae 
#' 
#' Reads a formula and derives pertinent information from iData object.
#' 
#' @param formula Object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData Object of class \code{\link{iData}} containing data represented in the provided formula.
#' @param impute Impute NA values (not yet implemented)
#' 
#' @return Creates list including:
#' \item{y} {Character value for the response variable (representative of an iGroup object within iData)}
#' \item{X} {The design matrix for the formula.}
#' \item{iData} {An iData object containing only pertinent information.}
#' 
#' @author Zachary P. Christensen
#' 
#' @example 
#' z <- iFormula(wb ~ Age, mydata)
#' out <- iModelMake(X = z$X, y = z$y, iData = z$iData)
#' out <- iModelSolve(out)
#' out <- summary(out, contrastMatrix = matrix(c(0, 1), nrow = 1), cthresh = 0)
#' 
#' @export iFormula
iFormula <- function(formula, iData, impute) {
  lhs <- formula[[2]]
  groups <- all.vars(lhs)
  for (i in seq_len(length(groups))) {
    if (!any(names(iData) == groups[i]))
      stop(paste(groups[i], " is not an iGroup found within iData", sep = ""))
  }
  
  # are there any NA values in demog 
  vars <- all.vars(formula[[3]])
  vartest <- TRUE
  for (i in seq_len(length(vars))) {
    vartest <- !any(is.na(iData@demog[vars[i]]))
    if (!vartest)
      break
  }
  
  # are there any images not common between groups
  grouptest <- TRUE
  tmpindex <- iData@index[groups]
  for (i in seq_len(nrow(tmpindex))) {
    grouptest <- !any(tmpindex[i, ] == 0)
    if (!grouptest)
      break
  }
  
  if (!vartest | !grouptest) # slower but gets rid of missing values
    iData <- select(iData, groups, vars, na.omit = TRUE)
  else # quick because should be using same pointers as original iData object
    iData <- select(iData, groups, vars, na.omit = FALSE)
  
  # probably not the most robust way to isolate right hand side
  rhs <- as.formula(paste("~", as.character(formula)[3], sep = ""))
  # This method drops out nonvariable terms (i.e. "-1")
  # tt <- terms(formula)
  # tt <- delete.response(tt)
  # rhs <- reformulate(attr(tt, "term.labels"))
  X <- model.matrix(rhs, data = iData@demog)
  list(y = groups, X = X, iData = iData)
}

#' Fit images to linear model
#' 
#' 
#' 
#' @param formula Object of class formula (or one that can be coerced to that class): a symbolic description of the modeto be fitted.
#' @param iData Object of class iData containing the variabels in the model.
#' @param weights An optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param optim Optimization method ("none", "REML", "IWLS").
#' @param its Maximum iterations for optimizing fitted models.
#' @param control A list of parameters for controlling the fitting process. See \code{\link{iControl}} for details.
#' @param verbose Enables verbose output. (default = \code{TRUE}).
#' 
#' @return \code{ilm} returns an object of class \code{\link{iModel}}.
#' 
#' @references 
#' Ashburner, Friston, and Penny (2004) Human Brain Function 2nd Edition.
#' 
#' @author Zachary P. Christensen
#' 
#' @example 
#' 
#' 
#' @export ilm
ilm <- function(formula, iData, weights = NULL, optim = "none", its, control, verbose = TRUE) {
  cl <- match.call()
  z <- iFormula(formula, iData)
  
  out <- list()
  for (i in seq_len(length(z$y))) {
    if (verbose)
      cat(paste("Fitting response", z$y[i], "\n"))
    out[[i]] <- iModelMake(X = z$X, y = z$y[i], iData = z$iData, weights = weights, control = control)
    out[[i]]@method <- c("ilm", optim)
    if (optim == "REML") {
      if (missing(its))
        its <- 32
      out[[i]]@xVi <- .estNonSphericity(object, its)
      out[[i]] <- iModelUpdate(out[[i]])
      out[[i]] <- iModelSolve(out[[i]])
    } else if (optim == "IWLS") {
      # compute leverages
      if (missing(its))
        its <- 200
      H <- diag(out[[i]]@X$X %*% tcrossprod(MASS::ginv(crossprod(out[[i]]@X$X), out[[i]]@X$X)))
      ores <- 1
      nres <- 10
      n <- 0
      while(max(abs(ores - nres)) > sqrt(1e-8)) {
        ores <- nres
        n <- n + 1
        
        if (n == 1) {
          W <- matrix(1, out[[i]]@dims$nimg, out[[i]]@dims$nvox)
          W[is.na((out[[i]]@iData@iList[[y]]@iMatrix[]))] <- 0
        }
        
        for (i in seq_len(out[[i]]@dims$nvox))
          out[[i]]@beta[, i] <- MASS::ginv(crossprod(out[[i]]@X$X, diag(W[, i])) %*% out[[i]]@X$X) %*% crossprod(out[[i]]@X$X, diag(W[, i])) %*% out[[i]]@iData@iList[[y]]@iMatrix[, i]
        
        if (n > its) {
          warning("ilm could not converge. Maximal number of iterations exceeded.");
          break
        }
        
        restmp <- out[[i]]@iData@iList[[y]]@iMatrix[] - out[[i]]@X$X %*% out[[i]]@beta[]
        
        mad <- rowMeans(abs(t(restmp) - colMeans(restmp)))
        restmp <- t(t(restmp) * (1 / mad))
        
        restmp <- restmp %*% H
        restmp <- abs(restmp) - control$os
        restmp(restmp < 0) <- 0
        nres <- sum(restmp(!is.na(restmp))^2)
        W <- (abs(restmp) < 1) %*% ((1 - restmp^2)^2)
        W(is.na(out[[i]]@iData@iList[[y]]@iMatrix[])) <- 0
        W(out[[i]]@iData@iList[[y]]@iMatrix[] == 0) <- 0
      }
      
      out[[i]]@res[] <- out[[i]]@iData@iList[[y]]@iMatrix[] - out[[i]]@X$X %*% out[[i]]@beta[]
      out[[i]]@mrss[1, ] <- colSums(out[[i]]@res[]^2) / out[[i]]@X$trRV
      
      if (verbose)
        cat(paste("Iterative reweighting finished after ", n, " iterations.", sep = ""))
    } else
      out[[i]] <- iModelSolve(out[[i]])
    
    if (out[[i]]@control$rft) {
      if (verbose)
        cat("Estimating FWHM/Resels. \n")
      smooth <- estSmooth(out[[i]]@res[], out[[i]]@iData@iList[[out[[i]]@y]]@mask, out[[i]]@X$rdf, scaleResid = FALSE, sample = out[[i]]@control$sar, verbose = verbose)
      out[[i]]@dims$fwhm <- smooth$fwhm
      out[[i]]@dims$rpvImage <- smooth[[2]]
      out[[i]]@dims$resels <- resels(out[[i]]@iData@iList[[out[[i]]@y]]@mask, smooth$fwhm)
    }
  }
  if (length(out) == 1)
    return(out[[1]])
  else
    return(structure(out, class = "multiModel"))
}

#' iModel Contrasts
#' 
#' @param object Object of class iModel.
#' @param contrastMatrix
#' @param cthresh Minimum desired cluster size (default = \code{150})
#' @param threshType A numeric value to threshold the statistical field or a character of the following methods:
#' \itemize{
#'	\item{cRFT} {Computes a threshold per expected cluster level probability.}
#'	\item{pRFT} {Uses the mask and pval calculates the minimum statistical threshold.}
#'	\item{cFDR} {Uses an uncorrected threshold at the alpha level and then computes and FDR threshold based on cluster maxima.}
#'	\item{pFDR} {Computes the fdr threshold for the entire field of voxels.}
#' }
#' @param pval The p-value for estimating the threshold (default = \code{0.05}).
#' @param pp The primary (initial) p-value for thresholding (only used for FDR methods; default = \code{0.001}).
#' @param verbose Enables verbose output.
#' 
#' @author Zachary P. Christensen
#' 
#' @name iModel-contrasts
NULL

#' @export
#' @docType methods
#' @description \strong{summary} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-contrasts
setMethod("anova", "iModel", function(object, contrastMatrix, cthresh = 150, threshType = "cFDR", pval = 0.05, pp = 0.001, verbose = TRUE) {
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(rep(1, ncol(object@X$X)), nrow = 1)
    colnames(contrastMatrix) <- colnames(object@X$X)
    if (any(colnames(contrastMatrix) == "(Intercept)"))
      contrastMatrix[, "Intercept"] <- 0
  }
  
  # check arguments----
  if (ncol(contrastMatrix) != object@dims$npred)
    stop("The contrast length must be equal to the number of predictors")
  ncon <- nrow(contrastMatrix)
  if (length(cthresh) != ncon)
    cthresh <- rep(cthresh[1], ncon)
  if (length(threshType) != ncon)
    threshType <- rep(threshType[1], ncon)
  if (length(pval) != ncon)
    pval <- rep(pval[1], ncon)
  if (length(pp) != ncon)
    pp <- rep(pp[1], ncon)
  
  
  # set isotropic option
  if (!object@control$iso)
    isotropic <- object@dims$rpvImage
  else
    isotropic <- NULL
  
  # set contrast names
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  
  
  mask <- object@iData@iList[[object@y]]@mask
  
  # ensure that if previous contrasts were set they won't be erased
  old <- length(object@C)
  conseq <- seq_len(ncon) + old
  
  ll <- object@control
  colnames(contrastMatrix) <- colnames(object@X$X)
  for (i in conseq) {
    c <- matrix(contrastMatrix[i - old, ])  # by default makes a 1 column matrix
    object@C[[i]] <- .setcon(connames[i - old], fieldType = "F", "c+", c, object@X$KWX)
    
    X1 <- .X1(object@C[[i]], object@X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)
    
    # calculate F-contrast----
    h <- .hsqr(object@C[[i]], object@X$KWX)
    ss <- (colSums((h %*% object@beta[])^2) / object@C[[i]]$dims$trMV)
    mrss <- object@mrss[]
    mrss[mrss == 0] <- 1
    object@C[[i]]$contrastImage <- makeImage(mask, ss / mrss)
    object@C[[i]]$cthresh <- cthresh[i - old]
    object@C[[i]]$threshType <- threshType[i - old]
    object@C[[i]]$pval <- pval[i - old]
    object@C[[i]]$pp <- pp[i - old]
    # summary of information----
    if (object@control$rft) {
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm,
                      c(object@C[[i]]$dims$idf, object@X$rdf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = object@C[[i]]$cthresh, object@C[[i]]$threshType, object@C[[i]]$pval, object@C[[i]]$pp, ll$n, NULL, verbose = verbose)
      if (z[1] == "NA") {
        object@C[[i]]$results <- "No voxels survive threshold."
        object@C[[i]]$sthresh <- 0
      } else {
        object@C[[i]]$sthresh <- z$threshold
        object@C[[i]]$clusterImage <- z$clusterImage
        object@C[[i]]$results$setLevel <- z$setLevel
        object@C[[i]]$results$peakLevel <- z$peakLevel
        object@C[[i]]$results$clusterLevel <- z$clusterLevel
      }
    } else {
      # don't use random field theory for summary----
      object@C[[i]]$sthresh <-
        statFieldThresh(object@C[[i]]$contrastImage, object@C[[i]]$pval, object@dims$nvox,
                        ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)),
                        c(object@C[[i]]$dims$idf, object$dims$rdf),
                        object@C[[i]]$fieldType, object@C[[i]]$threshType, object@C[[i]]$pp)
      
      object@C[[i]]$clusterImage <-
        labelClusters(object@C[[i]]$contrastImage, minClusterSize = object@C[[i]]$cthresh,
                      minThresh = object@C[[i]]$sthresh, maxThresh = Inf)
      
      object@C[[i]]$results <- labelStats(object@C[[i]]$contrastImage, object@C[[i]]$clusterImage)
    }
  }
  names(object@C)[conseq] <- connames
  return(object)
})

#' @export
#' @docType methods
#' @description \strong{summary} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-contrasts
setMethod("summary", "iModel", function(object, contrastMatrix, cthresh = 150, threshType = "pRFT", pval = 0.05, pp = 0.001, verbose = TRUE) {
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(rep(1, ncol(object@X$X)), nrow = 1)
    colnames(contrastMatrix) <- colnames(object@X$X)
    if (any(colnames(contrastMatrix) == "(Intercept)"))
      contrastMatrix[, "Intercept"] <- 0
  }
  
  # check arguments----
  if (ncol(contrastMatrix) != object@dims$npred)
    stop("The contrast length must be equal to the number of predictors")
  ncon <- nrow(contrastMatrix)
  if (length(cthresh) != ncon)
    cthresh <- rep(cthresh[1], ncon)
  if (length(threshType) != ncon)
    threshType <- rep(threshType[1], ncon)
  if (length(pval) != ncon)
    pval <- rep(pval[1], ncon)
  if (length(pp) != ncon)
    pp <- rep(pp[1], ncon)
  
  
  # set isotropic option
  if (!object@control$iso)
    isotropic <- object@dims$rpvImage
  else
    isotropic <- NULL
  
  # set contrast names
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  
  
  mask <- object@iData@iList[[object@y]]@mask
  
  # ensure that if previous contrasts were set they won't be erased
  old <- length(object@C)
  conseq <- seq_len(ncon) + old
  
  ll <- object@control
  colnames(contrastMatrix) <- colnames(object@X$X)
  for (i in conseq) {
    c <- matrix(contrastMatrix[i - old, ])  # by default makes a 1 column matrix
    object@C[[i]] <- .setcon(connames[i - old], fieldType = "T", "c+", c, object@X$KWX)
    
    X1 <- .X1(object@C[[i]], object@X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)
    
    # calculate t-contrast----
    Vc <- crossprod(object@C[[i]]$c, object@X$betaCov) %*% object@C[[i]]$c
    se <- sqrt(object@mrss[] * as.numeric(Vc))
    se[se == 0] <- 1  # prevent formation of NaN
    tvec <- crossprod(object@C[[i]]$c, object@beta[]) / se
    object@C[[i]]$contrastImage <- makeImage(mask, tvec)
    object@C[[i]]$cthresh <- cthresh[i - old]
    object@C[[i]]$threshType <- threshType[i - old]
    object@C[[i]]$pval <- pval[i - old]
    object@C[[i]]$pp <- pp[i - old]
    
    # summary of information----
    if (object@control$rft) {
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm,
                      c(object@C[[i]]$dims$idf, object@X$rdf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = object@C[[i]]$cthresh, object@C[[i]]$threshType, object@C[[i]]$pval, object@C[[i]]$pp, ll$n, NULL, verbose = verbose)
      if (z[1] == "NA") {
        object@C[[i]]$results <- "No voxels survive threshold."
        object@C[[i]]$sthresh <- 0
      } else {
        object@C[[i]]$sthresh <- z$threshold
        object@C[[i]]$clusterImage <- z$clusterImage
        object@C[[i]]$results$setLevel <- z$setLevel
        object@C[[i]]$results$peakLevel <- z$peakLevel
        object@C[[i]]$results$clusterLevel <- z$clusterLevel
      }
    } else {
      # don't use random field theory for summary----
      object@C[[i]]$sthresh <-
        statFieldThresh(object@C[[i]]$contrastImage, object@C[[i]]$pval, object@dims$nvox,
                        ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)),
                        c(object@C[[i]]$dims$idf, object$dims$rdf),
                        object@C[[i]]$fieldType, object@C[[i]]$threshType, object@C[[i]]$pp)
      
      object@C[[i]]$clusterImage <-
        labelClusters(object@C[[i]]$contrastImage, minClusterSize = object@C[[i]]$cthresh,
                      minThresh = object@C[[i]]$sthresh, maxThresh = Inf)
      
      object@C[[i]]$results <- labelStats(object@C[[i]]$contrastImage, object@C[[i]]$clusterImage)
    }
  }
  names(object@C)[conseq] <- connames
  return(object)
})