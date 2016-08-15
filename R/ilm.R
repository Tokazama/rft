# to do:
# multiple modality eigenanatomy stuff

#' iData Model Formulae 
#' 
#' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param iData an object of class \code{\link{iData}} containing data represented in the provided formula
#' @param ...
#' 
#' 
#' 
#' @export iFormula
iFormula <- function(formula, iData, ...) {
  lhs <- formula[[2]]
  groups <- all.vars(lhs)
  for (i in seq_len(length(groups))) {
    if (!any(names(iData) == groups[i]))
      stop(paste(groups[i], " is not an iGroup found within iData", sep = ""))
  }
  cleanData <- na.omit(iData, groups, all.vars(formula[[3]]))
  
  tt <- terms(formula)
  tt <- delete.response(tt)
  rhs <- reformulate(attr(tt, "term.labels"))
  X <- model.matrix(rhs, data = cleanData@demog)
  list(y = groups, X = X)
}

#' Fit images to linear model
#' 
#' 
#' 
#' @param formula Object of class formula (or one that can be coerced to that class): a symbolic description of the modeto be fitted.
#' @param iData Object of class iData containing the variabels in the model.
#' @param weights 
#' @param optim optimization method ("none", "REML", "IWLS")
#' @param control 
#' @param verbose 
#'
#'
#' @export ilm
ilm <- function(formula, iData, weights = NULL, optim = "none", control, verbose = TRUE) {
  cl <- match.call()
  z <- iFormula(formula, iData)
  
  ctrl <- iControl()
  if (!missing(control))
    ctrl[names(control)] <- control
  
  out <- list()
  for (i in seq_len(length(z$y))) {
    if (verbose)
      cat(paste("Fitting response", z$y[i], "\n"))
    out[[i]] <- iModel(X = z$X, y = z$y[i], data = iData, weights = weights, control = control, ...)
    out[[i]]@method <- c("ilm", optim)
    if (optim == "REML") {
      out[[i]]@xVi <- estNonSphericity(object)
      out[[i]] <- iModelUpdate(out[[i]])
      out[[i]] <- iModelSolve(out[[i]])
    } else if (optim == "IWLS") {
      # compute leverages
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
        
        if (n > out[[i]]@control$mi) {
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
    
    if (ctrl$rft) {
      smooth <- estSmooth(out[[i]]@res[], out[[i]]@iData[[out[[i]]@y]]@mask, out[[i]]@X$rdf, scaleResid = FALSE, sample = control$sar, verbose = verbose)
      out[[i]]@dims$fwhm <- smooth$fwhm
      out[[i]]@dims$rpvImage <- smooth$rpvImage
      out[[i]]@dims$resels <- resels(out[[i]]@iData[[x@y]]@mask, smooth$fwhm)
    }
  }
  if (length(out) == 1)
    return(out[[1]])
  else
    return(structure(out, class = "multiModel"))
}