# to do:
# multiple modality eigenanatomy stuff
# impute for ilm

#' iData Model Formulae 
#' 
#' Reads a formula and derives pertinent information from iData object. If any 
#' variables in the formula have \code{NA} values a new \code{\link{iData}}
#' object is generated either omitting those values or if function is provided
#' for \code{impute} argument \code{NA} values are imputed (see 
#' \code{\link{antsrimpute}}). Unlike \code{\link{antsrimpute}} there's no
#' options for additional arguments to be passed to the impute function.
#' Therefore, if additional arguments are necessary one may need to create a
#' temporary function including those arguments as default (see examples).
#' 
#' @param formula Object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData Object of class \code{\link{iData}} containing data represented in the provided formula.
#' @param impute Function for imputing NA values.
#' 
#' @return Creates list including:
#' \item{y} {Character value for the response variable (representative of an iGroup object within iData)}
#' \item{X} {The design matrix for the formula.}
#' \item{iData} {An iData object containing only pertinent information.}
#' 
#' @author Zachary P. Christensen
#' 
#' @example 
#' 
#' ilist <- getANTsRData("population")
#' mask <- getMask(ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup1 <- iGroup(imat, "pop1", mask, modality = "T1")
#' 
#' 
#' ilist <- lappend(ilist, ilist[[1]])
#' imat <- imageListToMatrix(ilist, mask)
#' iGroup2 <- iGroup(imat, "pop2", mask, modality = "T1")
#' 
#' demog <- data.frame(id = c("A", "B", "C", NA),
#'   age = c(11, 7, 18, 22), sex = c("M", "M", "F", "F"))
#'   
#' bool1 <- c(TRUE, TRUE, TRUE, FALSE)
#' bool2 <- c(TRUE, TRUE, TRUE, TRUE)
#' 
#' # create iData object that holds demographics info
#' mydata <- iData(list(iGroup1, iGroup2), c(bool1, bool2), demog)
#' 
#' z <- iFormula(iGroup1 ~ age, mydata)
#' 
#' 
#' # quick function for mean with custom defaults
#' myfunc <- function(x) {
#'   mean(x, trim = .1)
#' }
#' 
#' z <- iFormula(iGroup1 ~ age, mydata, myfunc)
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
  
  if (!vartest && !missing(impute)) {
    iData@demog <- antsrimpute(iData@demog[vars], impute)
    vartest <- TRUE
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
      out[[i]]@xVi <- .estNonSphericity(out[[i]], its)
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
      W <- matrix(1, out[[i]]@dims$nimg, out[[i]]@dims$nvox)
      while(max(abs(ores - nres)) > (1e-4)) {
        ores <- nres
        n <- n + 1
        
        if (n > its) {
          warning("ilm could not converge. Maximal number of iterations exceeded.");
          break
        }
        
        nres <- 0
        for (j in seq_len(out[[i]]@dims$nvox)) {
          if (i == 1)
            W[, j][is.na(out[[i]]@iData@iList[[out[[i]]@y]]@iMatrix[, j])] <- 0
          
          out[[i]]@beta[, j] <- MASS::ginv(crossprod(out[[i]]@X$X, diag(W[, j])) %*% out[[i]]@X$X) %*% crossprod(out[[i]]@X$X, diag(W[, j])) %*% out[[i]]@iData@iList[[out[[i]]@y]]@iMatrix[, j]
          out[[i]]@res[, j] <- out[[i]]@iData@iList[[out[[i]]@y]]@iMatrix[, j] - out[[i]]@X$X %*% out[[i]]@beta[, j]
          restmp <- out[[i]]@res[, j]
          mad <- mean(abs(restmp - mean(restmp)))
          restmp <- restmp / mad
          
          # need to figure this one out
          restmp <- restmp * H
          
          restmp <- abs(restmp) - control$os
          restmp[restmp < 0] <- 0
          
          wtmp <- (abs(restmp) < 1) * ((1 - restmp^2)^2)
          wtmp[is.na(out[[i]]@iData@iList[[out[[i]]@y]]@iMatrix[, j])] <- 0
          wtmp[out[[i]]@iData@iList[[out[[i]]@y]]@iMatrix[, j] == 0] <- 0
          W[, j] <- wtmp
          nres <- nres + sum(restmp[!is.na(restmp)]^2)
        }
      }
      
      out[[i]]@res[] <- out[[i]]@iData@iList[[out[[i]]@y]]@iMatrix[] - out[[i]]@X$X %*% out[[i]]@beta[]
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
  if (any(is.na(contrastMatrix)))
    stop("There cannot be any NA values in the contrastMatrix.")
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
    if (verbose)
      cat("Calculating contrast: ", connames[i - old], ".\n")
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
                      k = object@C[[i]]$cthresh, object@C[[i]]$threshType, object@C[[i]]$pval, object@C[[i]]$pp, ll$n, verbose = verbose)
      if (is.null(z)) {
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
  if (any(is.na(contrastMatrix)))
    stop("There cannot be any NA values in the contrastMatrix.")
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
    if (verbose)
      cat("Calculating contrast: ", connames[i - old], ".\n")
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
                      k = object@C[[i]]$cthresh, object@C[[i]]$threshType, object@C[[i]]$pval, object@C[[i]]$pp, ll$n, verbose = verbose)
      if (is.null(z)) {
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