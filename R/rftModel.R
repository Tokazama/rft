#' Create a random field object
#'
#' @param X design matrix [subject (n) x predictor (p)]
#' @param y response matrix (typically an image matrix) [n x voxels (v)]
#' @param mask a mask of class antsImage
#' @param contrastMatrix matrix specifying contrasts [contrasts x p]
#' @param weights vector of weights corresponding to each row of the design matrix and response matrix or 
#' @param method the method used to fit the model (default = "OLS")
#' @param findVar logical. calculate the response matrix covariance
#' @param findResels logical. estimate the resolution elements of the fitted model
#' @param control a list of parameters for controlling the fitting process. See rftControl documentation for more details.
#' @param verbose activates verbose output
#' 
#' @return
#' \item{call} {the matched call to the original argument}
#' \item{modParams} {model parameters}
#' \itemize{
#'   \item{X} {design matrix}
#'   \item{R} {the residual projection matrix}
#'   }
#' \item{varParams} {variance parameters}
#' \itemize{
#'   \item{V} {}
#'   \item{Cy} {}
#'   }
#' \item{rss} {residual sum of squares}
#' \item{coefficients} {a matrix of predictors by voxels}
#' \item{dof} {degrees of freedom [effective degreesof freedom, residual degrees of freedom]}
#' \item{dims} {terms describing the dimensions of the fitted model}
#' \itemize{
#'   \item{rpvImage} {resel per voxel image}
#'   \item{fwhm} {full-widths at half-maxima}
#'   \item{resels} {resolution elements}
#'   \item{nv} {number of voxels analyzed}
#'   \item{ns} {number of subjects}
#'   \item{np} {number of predictors}
#'   }
#' \item{contrastMatrix} {the input matrix argument with labels}
#' \item{control} {a list of the control parameters}
#' \item{statdir} {the directory where results are saved}
#' \item{method} {the method used to fit the model}
#' 
#' 
#' @export rftModel
rftModel <- function(X, y, mask, contrastMatrix, weights = diag(nrow(X)),
                     method = "OLS", findVar = FALSE, findResels = FALSE,
                     control = list(), verbose = FALSE) {
  # check parameters-----------------------------------------------------------
  n  <- nrow(X)
  p  <- ncol(X)
  v  <- ncol(y)
  controlvals <- rftControl()
  if (!missing(control))
    controlvals[names(control)] <- control
  
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(1, 1, p)
  } else {
    if (ncol(contrastMatrix) != p)
      stop("Number of columns in the design matrix must equal the number of columns in contrastMatrix.")
  }
  wdim <- dims(weights)
  if (is.null(wdim)) {
    weights <- diag(weights)
    wdim <- dim(weights)
  }
  if (wdim[1] != n)
    stop("the number of weights must equal the sample size")
  
  # satterthwhite approximation------------------------------------------------
  # R <- diag(n) - X %*% MASS::ginv(crossprod(X, V) %*% X) %*% crossprod(X, V)
  R <- diag(n) - X %*% MASS::ginv(crossprod(X, weights), %*% crossprod(X, weights))
  RV <- R %*% weights
  trRV <- sum(diag(RV))
  trRVRV <- sum(diag(RV %*% RV))
  
  c <- matrix(1, 1, p)
  C0 <- diag(p) - c %*% MASS::ginv(c)
  X0 <- X %*% C0
  R0 <- diag(n) - X0 %*% MASS::ginv(X0)
  M  <- R0 - R
  MV <- M %*% weights
  
  dof <- c((trMV^2) / sum(diag(MV %*% MV)), (trRV^2) / sum(diag(RV %*% RV)))
  # OLS fit--------------------------------------------------------------------
  if (verbose)
    cat("Fitting model")
  XV <- X %*% weights
  XX <- MASS:ginv(crossprod(XV, X))
  YV <- y %*% weights
  beta <- XX %*% crossprod(XV, YV)
  
  fv <- beta %*% X
  r <- y - fv
  
  fthresh <- qf(1 - controlvals$criticalF, dof[1], dof[2])
  rss <- colSums(r^2)
  fstat <- (((fv - colMeans(y))^2) / trMV) / (rss / trRV)
  good <- cov(y[fstat > fthresh])
  
  # pooled variance------------------------------------------------------------
  if (findVar) {
    if (verbose)
      cat("Calculating pooled variance \n")
    y <- y[, good] %*% diag(trMV / rss[good])
    Cy <- tcrossprod(y)
  }
  
  # estimate fwhm/resel--------------------------------------------------------
  if (findResels == TRUE) {
    if (verbose)
      cat("Estimating FWHM. \n")
    mysmooth <- estSmooth(r, newmask, rdf, scaleResid = TRUE, sample = sample, verbose = verbose)
    if (verbose)
      cat("Estimating resels. \n")
    myresels <- resels(newmask, mysmooth$fwhm)
  } else {
    mysmooth <- NULL
    myresels <- NULL
  }
  
  # create rftModel object-----------------------------------------------------
  structure(list(call = NULL,
                 modParams = list(X = X, R = R),
                 varParams = list(V = V, Cy = Cy),
                 rss = rss,
                 coefficients = beta,
                 dof = c(df1, df2),
                 dims = list(mask = newmask,
                             rpvImage = mysmooth$rpvImage,
                             fwhm = mysmooth$fwhm,
                             resels = myresels,
                             nv = v,
                             ns = n,
                             np = p),
                 contrastMatrix = contrastMatrix,
                 control = controlvals,
                 statdir = statdir,
                 method = method),
            class = "rftModel")
}


#' @export update.rftModel
update.rftModel <- function(object, X, y, mask, contrastMatrix, weights, findResels = FALSE, method) {
  if (!missing(X) | !missing(y) | !missing(weights) | findResels = TRUE) {
    call <- object$modParams$call
    extras <- match.call()
    if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing])
        call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    if (evaluate)
      eval(call, parent.frame())
    else 
      call
  } else {
    if (!missing(mask)) {
      object$dims$mask <- antsImageClone(mask)
    }
    if (!missing(contrastMatrix)) {
      object$contrastMatrix <- contrastMatrix
    }
    object
  }
}

#' @export print.rftModel
print.rftModel <- function(x) {
  cat("Call: \n")
  cat(x$call, "\n\n")
  
  cat("Coefficients: \n")
  ncoef <- nrow(x$beta)
  for (i in 1:ncoef) {
    cat("  ", colnames(x$coefficients)[i], "\n")
  }
  cat("\n")
  
  cat("Degrees of freedom = ", x$dof, "\n")
  cat("Voxels = ", x$dims$vn, "\n")
  cat("FWHM = ", x$dims$fwhm, "\n")
  cat("Resels = ", x$dims$resels, "\n\n")
  
  cat("Directory: ", x$statdir, "\n")
}

#' @export residuals.rftModel
residuals.rftModel <- function(object, scaleResid = TRUE) {
  call <- object$modParams$call
  z <- do.call(".lm.fit", list(x = object$modParams$X, y = call$y), parent.frame())$residuals
  if (scaleResid == TRUE)
    z / rss
  else
    z
}
