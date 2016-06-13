#' Create an object of class rftModel
#'
#' @param X design matrix [subject (n) x predictor (p)]
#' @param y n by p input image matrix, where n is the number of subjects and p is the number of voxels.
#' @param mask Mask image of type antsImage.
#' @param contrastMatrix matrix specifying contrasts [contrasts x p]
#' @param weights vector of weights corresponding to each row of the design matrix and response matrix
#' @param V error covariance matrix (defaults to diag(n))
#' @param method the method used to fit the model (default = "OLS")
#' @param findVar logical. calculate the response matrix covariance
#' @param findResels logical. estimate the resolution elements of the fitted model
#' @param control a list of parameters for controlling the fitting process. See rftControl documentation for more details.
#' @param verbose activates verbose output
#' 
#' @return list(X = X, XX = pWX, trRV = trRV, df2 = df2)
#' \item{call} {the matched call to the original argument}
#' \item{modParams} {model parameters}
#' \itemize{
#'   \item{X} {design matrix}
#'   \item{XX} {}
#'   \item{trRV} {}
#'   \item{df2} {residual degrees of freedom for the model}
#'   }
#' \item{varParams} {variance parameters}
#' \itemize{
#'   \item{V} {}
#'   \item{Cy} {error covariance matrix}
#'   }
#' \item{residuals}
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
#' \item{contrasts} {list of contrast parameters}
#' \itemize{
#'   \item{contrast} {the contrast}
#'   \item{trMV}
#'   \item{df1} {the effective degree of freedom}
#'   }
#' \item{control} {list of the control parameters}
#' \item{statdir} {the directory where results are saved}
#' \item{method} {the method used to fit the model}
#' 
#' contrast = c, trMV = trMV, df1 = df1
#' @export rftModel
rftModel <- function(X, y, mask, contrastMatrix = NULL, weights = NULL, V = NULL, control = list(), statdir = NULL, verbose = TRUE) {
  # check parameters-----------------------------------------------------------
  n  <- nrow(X)
  p  <- ncol(X)
  v  <- ncol(y)
  cl <- match.call()
  controlvals <- rftControl()
  if (!missing(control))
    controlvals[names(control)] <- control
  
  # set contrasts--------------------------------------------------------------
  predNames <- colnames(X)
  modcon <- matrix(1, 1, p)
  if (any(predNames == "Intercept"))
    modcon[1,1] <- 0
  if (is.null(contrastMatrix)) {
    contrastMatrix <- modcon
    conNames <- "model_contrast"
  } else {
    condims <- dim(contrastMatrix)
    if (is.null(condims)) {
      contrastMatrix <- matrix(contrastMatrix, 1, p)
      condims <- dim(contrastMatrix)
    }
    if (condims[2] != p)
      stop("Number of columns in the design matrix must equal the number of columns in contrastMatrix.")
    conNames <- rownames(contrastMatrix)
    if (is.null(conNames))
      conNames <- paste("contrast_", 1:condims[1], sep = "")
    conNames <- c("model_contrast", conNames)
    contrastMatrix <- rbind(modcon, contrastMatrix)
  }
  ncon <- nrow(contrastMatrix)
  colnames(contrastMatrix) <- predNames
  rownames(contrastMatrix) <- conNames
  
  # set weights/error-variance parameters--------------------------------------
  if (is.null(V))
    V <- diag(n)
  if (is.null(weights)) {
    iV <- sqrt(MASS::ginv(V))
    weights <- iV * (abs(iV) > 1e-6)
  } else {
    if (!missing(weights)) {
      wdim <- dims(weights)
      if (is.null(wdim)) {
        weights <- diag(weights)
        wdim <- dim(weights)
      }
      if (wdim[1] != n)
        stop("The number of weights must equal the sample size")
    }
  }
  
  # account for weights in variance and design matrix
  WX <- weights %*% X  # weighted design matrix
  pWX <- MASS::ginv(crossprod(WX))  # weighted design matrix projector
  V <- weights %*% tcrossprod(V, weights)
  
  # degrees of freedom---------------------------------------------------------
  u <- svd(WX)$u
  Vu <- V %*% u
  trV <- sum(diag(V))
  trRVRV <- sum(diag(crossprod(V)))
  trRVRV <- trRVRV - 2 * sum(diag(crossprod(Vu)))
  trRVRV <- trRVRV + sum(diag(crossprod(crossprod(u, Vu))))
  trMV <- sum(u * Vu)
  trRV <- trV - trMV
  df2 <- (trRV^2) / trRVRV
  
  contrasts <- list()
  for (i in 1:ncon) {
    c <- matrix(contrastMatrix[i, ], p, 1)
    C0 <- diag(p) - c %*% MASS::ginv(c)
    X0 <- X %*% C0
    u <- svd(X0)$u
    Vu <- V %*% u
    trMV <- sum(u * Vu)
    trMVMV <- sum(diag(crossprod(crossprod(u, Vu))))
    df1 <- (trMV^2) / trMVMV
    contrasts[[i]] <- list(contrast = c, trMV = trMV, df1 = df1)
  }
  names(contrasts) <- conNames
  
  # OLS fit--------------------------------------------------------------------
  if (verbose)
    cat("Fitting model \n")
  WY <- weights %*% y
  olsfit <- .lm.fit(WX, WY)
  rss <- colSums(olsfit$residuals^2)
  
  # pooled variance------------------------------------------------------------
  if (verbose)
    cat("Calculating pooled variance \n")
  
  fthresh <- qf(1 - controlvals$criticalF, contrasts[[1]]$df1, df2)
  se <- sqrt((rss / trRV) * (crossprod(contrasts[[1]]$contrast, pWX) %*% contrasts[[1]]$contrast))
  se[se == 0] <- 1
  fstat <- (crossprod(contrasts[[1]]$contrast, olsfit$coefficients) / se)^2
  
  #good <- fstat > 0
  good <- fstat > fthresh
  
  tmpy <- y[, good]
  q <- sqrt(contrasts[[1]]$trMV / rss[good])
  ntmpy <- ncol(tmpy)
  
  if (verbose)
    progress <- txtProgressBar(min = 0, max = i, style = 3)
  for (i in 1:ntmpy) {
    tmpy[, i] <- tmpy[, i] * q[i]
    if (verbose)
      setTxtProgressBar(progress, i)
  }
  if (verbose)
    close(progress)
  
  Cy <- tcrossprod(tmpy)
  
  rownames(olsfit$coefficients) <- predNames
  
  # create rftModel object-----------------------------------------------------
  structure(list(rftModelCall = cl,
                 modParams = list(X = X, XX = pWX, trRV = trRV, df2 = df2),
                 varParams = list(V = V, Cy = Cy),
                 residuals = olsfit$residuals,
                 rss = rss,
                 coefficients = olsfit$coefficients,
                 dims = list(mask = antsImageClone(mask),
                             rpvImage = NULL,
                             fwhm = NULL,
                             resels = NULL,
                             nv = v,
                             ns = n,
                             np = p),
                 contrasts = contrasts,
                 control = controlvals,
                 statdir = statdir),
            class = "rftModel")
}


#' @export update.rftModel
rftModelUpdate <- function(object, X, y, mask, contrastMatrix, weights, V, control, statdir, verbose) {
  if (!missing(X) | !missing(y) | !missing(weights) | !missing(V)) {
    call <- object$rftModelCall
    extras <- match.call()
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    object <- eval(call, parent.frame())
    if (verbose)
      cat("Estimating FWHM. \n")
    object$dims[c("fwhm", "rpvImage")] <- estSmooth(object$residuals, object$dims$mask, object$modParams$df2, scaleResid = object$control$scaleResid, sample = object$control$sampleResid, verbose = object$verbose)
    if (verbose)
      cat("Estimating resels. \n")
    object$dims$resels <- resels(object$dims$mask, object$dims$fwhm)
  } else {
    if (!missing(statdir))
      object$statdir <- statdir
    if (!verbose)
      object$verbose <- verbose
    if (!missing(control))
      object$control[names(control)] <- control
    if (!missing(mask))
      object$dims$mask <- antsImageClone(mask)
    
    # set contrasts--------------------------------------------------------------
    if (!missing(contrastMatrix)) {
      predNames <- colnames(X)
      modcon <- matrix(1, 1, p)
      if (any(predNames == "Intercept"))
        modcon[1,1] <- 0
      if (is.null(contrastMatrix)) {
        contrastMatrix <- modcon
        conNames <- "model_contrast"
      } else {
        condims <- dim(contrastMatrix)
        if (is.null(condims)) {
          contrastMatrix <- matrix(contrastMatrix, 1, p)
          condims <- dim(contrastMatrix)
        }
        if (condims[2] != p)
          stop("Number of columns in the design matrix must equal the number of columns in contrastMatrix.")
        conNames <- rownames(contrastMatrix)
        if (is.null(conNames))
          conNames <- paste("contrast_", 1:condims[1], sep = "")
        conNames <- c("model_contrast", conNames)
        contrastMatrix <- rbind(modcon, contrastMatrix)
      }
      ncon <- nrow(contrastMatrix)
      colnames(contrastMatrix) <- predNames
      rownames(contrastMatrix) <- conNames
      
      # account for weights in variance and design matrix
      WX <- weights %*% X  # weighted design matrix
      pWX <- MASS::ginv(crossprod(WX))  # weighted design matrix projector
      V <- weights %*% tcrossprod(V, weights)
      
      # degrees of freedom---------------------------------------------------------
      u <- svd(WX)$u
      Vu <- V %*% u
      trV <- sum(diag(V))
      trRVRV <- sum(diag(crossprod(V)))
      trRVRV <- trRVRV - 2 * sum(diag(crossprod(Vu)))
      trRVRV <- trRVRV + sum(diag(crossprod(crossprod(u, Vu))))
      trMV <- sum(u * Vu)
      trRV <- trV - trMV
      df2 <- (trRV^2) / trRVRV
      
      contrasts <- list()
      for (i in 1:ncon) {
        c <- matrix(contrastMatrix[i, ], p, 1)
        C0 <- diag(p) - c %*% MASS::ginv(c)
        X0 <- X %*% C0
        u <- svd(X0)$u
        Vu <- V %*% u
        trMV <- sum(u * Vu)
        trMVMV <- sum(diag(crossprod(crossprod(u, Vu))))
        df1 <- (trMV^2) / trMVMV
        contrasts[[i]] <- list(contrast = c, trMV = trMV, df1 = df1)
      }
      names(contrasts) <- conNames
      object$contrasts <- contrasts
    }
  }
  # update models with resels--------------------------------------------------
  if (is.null(object$dims$fwhm)) {
    if (verbose)
      cat("Estimating FWHM. \n")
    object$dims[c("fwhm", "rpvImage")] <- estSmooth(object$residuals, object$dims$mask, object$modParams$df2, scaleResid = object$control$scaleResid, sample = object$control$sampleResid, verbose = object$verbose)
    if (verbose)
      cat("Estimating resels. \n")
    object$dims$resels <- resels(object$dims$mask, object$dims$fwhm)
  }
  object
}

#' @export print.rftModel
print.rftModel <- function(x) {
  cat("Call: \n")
  print(x$rftModelCall)
  cat("\n")
  cat("Coefficients: \n")
  ncoef <- nrow(x$coefficients)
  for (i in 1:ncoef) {
    cat("  ", rownames(x$coefficients)[i], "\n")
  }
  cat("\n")
  
  cat("Effictive degrees of freedom = ", x$contrasts[[1]]$df1, "\n")
  cat("Residual degrees of freedom = ", x$modParams$df2, "\n")
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
