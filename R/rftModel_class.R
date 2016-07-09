## This file contains the following functions pertinant to objects of class rftModel (in order):
##
## coef              (done; check)
## fitted            (done; check)
## model.matrix      (done; check)
## residuals         (done; check)
## rftControl        (IN PROGRESS; check)
## rftModel          (IN PROGRESS; check)
## update            (done; check)
## weights           (done; check)

# TO DO:
#   use @method print rftModel
#   or
#   use @S3method print rftModel



##' Class "rftModel" of Random Field Theory Fitted Models
##'
##' A fitted random field theory model is represented as a rftModel object.
##'
##' @name rftModel-class
##'
##' @docType class
##'
##' @section Objects from the Class: Objects are created by calls to
##' \code{\link{rftLm}}
##'
##' @seealso
##'
##' @keywords classes
##' @examples
##'
##' showClass("rftModel")
##' methods(class="rftModel")
##' @export

#' @importFrom stats coef
#' @S3method coef rftModel
coef.rftModel <- function(object) {
  object$beta
}

## @method confint rftModel
## @describeIn rftModel-class
# setMethod("confint", "rftModel", function(object, parm, level = 0.95, ...) {
#   cf <- object$B
#   pnames <- rownames(cf)
#   if (missing(parm))
#     parm <- pnames
#   else if (is.numeric(parm))
#     parm <- pnames[parm]
#
#   a <- (1 - level)/2
#   a <- c(a, 1 - a)
#
#   pct <- format.perc(a, 3)
#
#   fac <- qnorm(a)
#
#   ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
#
#   ses <- sqrt(diag(vcov(object)))[parm]
#   ci[] <- cf[parm] + ses %o% fac
#   ci
#   })


## @method deviance rftModel
## @describeIn rftModel-class
## setMethod("deviance", "rftModel", function(object) {
#
# })

#' @importFrom stats fitted
#' @S3method fitted rftModel
fitted.rftModel <- function(object) {
  object$imgData[object$y]$imageMatrix - object$beta %*% object$X
}

#' @importFrom stats model.matrix
#' @S3method model.matrix rftModel
model.matrix.rftModel <- function(object) {
  object$X
}

#' @importFrom stats residuals
#' @S3method residuals rftModel
residuals.rftModel <- function(object) {
  object$resid
}

#' Control parameters for RFT based analyses
#'
#' @param maxIter
#' @param scaleResid
#' @param sampleResid
#' @param criticlF critical F-threshold for selecting voxels over which the non-sphericity is estimated (default = 0.001)
#' @param threshType
#' @param threshPval
#' @param fdrThreshPval
#' @param conjuncImages
#' @param isotropic
#'
#' @export rftModelControl
rftControl <- function() {
  list(criticalF = 0.05,
       maxIter = 25,
       scaleResid = TRUE,
       sampleResid = 64,
       chunksize = function(nimg) ((2^23) / nimg),
       epsilon = 1e-08)  # for the 'tol' argument in .lm.fit
}

#' Create an object of class rftModel
#'
#'
#'
#'
#'
#' @field imgData
#' @field y name of image group from imgData used in fitting rftModel
#' @field beta
#' @field resid
#' @field mrss
#' @field X
#' @field XV correlation after \code{W} is applied
#' @field XX
#' @field W
#' @field KWX
#' @field Vi
#'
#' @field K
#' @field rpvImage
#' @field dims dimensions describing the model:
#' \itemize{
#'  \item{nimg} {number of images}
#'  \item{npred} {number of predictors}
#'  \item{trRV} {trace of R %*% V}
#'  \item{trRVRV} {trace of RV %*% RV}
#'  \item{rdf} {residual degrees of freedom}
#' }
#' @field control
#' @field call
#'
#' @export rftModel-class
rftModel <-
  setRefClass("rftModel",
              fields = list(
                imgData = "modData",
                y = "character",
                X = "matrix",
                W = "matrix",
                formula = "formula",
                B = "matrix",
                CB = "matrix",
                R = "matrix",
                MRSS = "matrix",
                XV = "matrix",
                XX = "matrix",
                KWX = "list",
                xvi = "list",
                K = "ANY",
                CY = "matrix",
                dims = "list",
                d = "list",  # design information
                rpvImage = "antsImage",
                control = "list",
                call = "call"),
              methods = list(
                setB = function() {
                  'set the model coefficients/beta'
                  B <<- XX %*% rftFilter(K, diag(W) * imgData$imgList[y]$imageMatrix)
                },
                setCB = function() {
                  CB <<- XX %*% tcrossprod(XV, XX)
                },
                setCY = function(cf) {
                  'spatially whitened tcrossprod(y) (used by ReML to estimate hyperparameters)'
                  if (missing(cf))
                    cf <- 0.05
                  CY < util_cy(.self, cf)
                },
                setCall = function(call) {
                  'set call field'
                  call <<- call
                },
                setTrRV = function(oneout = FALSE) {
                  'traces for effective df calculation. If oneout = TRUE then returns trRV.'
                  rk <- KWX$rk
                  sL <- nrow(KWX$X)

                  u <- KWX$u[, seq_len(rk)]
                  if (oneout) {
                    if (rk == 0)
                      return(0)
                    else
                      trmv <- sum(u %*% (XV %*% u))
                    return(sum(diag(XV)) - trmv)
                  } else {
                    if (rk == 0) {
                      trmv <- 0
                      tmp <- norm(XV, "f")^2
                      trv <- sum(diag(XV))
                    } else {
                      Vu <- XV %*% u
                      trv <- sum(diag(XV))
                      tmp <- norm(XV, "F")^2
                      tmp <- tmp - 2 * norm(Vu, "f")^2
                      dims$trRVRV <<- tmp + norm(crossprod(u, Vu), "f")^2
                      trmv <- sum(u %*% Vu)
                    }
                    dims$trRV <<- trV - trmv
                  }
                  dims$rdf <<- (trRV^2) / trRVRV
                  dims$npred <<- ncol(X)
                  dims$idf <<- ncol(X) - 1
                },
                getResid = function() {
                  'return residual forming matrix or set residuals'
                  KWY <- rftFilter(K, diag(W) * imgData$imgList[y]$imageMatrix)

                  if (KWX$rk < nrow(KWX$X) - KWX$rk) {
                    out <- KWY - KWX$u[, seq_len(KWX$rk)] %*%
                      crossprod(KWX$u[, seq_len(KWX$rk)], KWY)
                  } else {
                    if (ncol(KW) < 5 * ncol(KWX$X)) {
                      R <- diag(ncol(KWX$X)) - crossprod(KWX$v[, seq_len(KWX$rk)], KWX$v[, seq_len(KWX$rk)])
                      out <- R %*% KWY
                    } else {
                      out  <- list()
                      out$X <- t(KWX$X)
                      out$v <- KWX$u
                      out$u <- KWX$v

                      out$oP <- KWX$opp
                      out$oPp <- KWX$op

                      dimX <- dim(out$X)
                      if (x$rk > 0) {
                        if (dimX[1] >= dimX[2]) {
                          if (out$rk == dimX[2])
                            out <- matrix(0, dimX[2], 1)
                          else
                            out <- out$v[, out$rk:dimX[2]]
                        } else
                          out <- MASS::Null(out$X)
                      } else
                        out <- diag(dimX[2])
                      out <- n %*% crossprod(n, KWY)
                    }
                  }
                  resid <<- out
                },
                setResels = function (verbose = NULL) {
                  'estimate FWHM and set resels'
                  tmp <- colSums(resid^2)
                  mrss <<- tmp / dims$trRV
                  smooth <- estSmooth(resid, modData$mask, dims$rdf, scaleResid = FALSE,
                                      control$sampleResid, verbose)
                  dims$fwhm <<- smooth$fwhm
                  rpvImage <<- smooth$rpvImage
                  dims$resels <<- resels(modData$mask, smooth$fwhm)
                },
                setW = function(weights) {
                  'set the weight matrix'
                  if (missing(weights)) {
                    iV <- sqrt(MASS::ginv(Vi$V))
                    weights <- iV * (abs(iV) > 1e-6)
                  } else {
                    if (class(weights) == "numeric" && length(weights) == nimg)
                      W <<- diag(weights)
                    else if (all(dim(weights) == nimg))
                      W <<- weights
                    else
                      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
                  }
                },
                setXX = function(x, y, check = FALSE) {
                  'set the pseudoinverse of KWX'
                  if (KWX$rk > 0 )
                    out <- KWX$v[, seq_len(KWX$rk)] %*%
                      tcrossprod(diag(rep(1, KWX$rk) / KWX$d[seq_len(KWX$rk)]),
                                 KWX$u[, seq_len(KWX$rk)])
                  else
                    out <- matrix(0, ncol(KWX$X), nrow(KWX$X))
                  if (!missing(y))
                    out <- out %*% y
                  if (check)
                    out[abs(out) < x$tol] <- 0
                  XX <<- out
                },
                setXV = function() {
                  'set correlations after weighting and filtering are applied'
                  XV <<- rftFilter(K, t(rftFilter(K, W %*% Vi$V %*% t(W))))
                },
                setKWX = function() {
                  'set the filtered and whitened design matrix'
                  tmp <- rftFilter(K, diag(W) * X)
                  out <- svd(tmp)
                  out$X <- tmp
                  out$tol <- max(dim(X)) * max(absx$d) * .Machine$double.eps
                  out$rk <- sum(x$d > x$tol)
                  KWX <<- out
                },
                show = function() {
                  'method for printing rftModel'
                  cat("Random Field Theory model fitted by ", call[[1]])
                  cat("Call: \n")
                  print(call)

                  cat("\nCoefficients: \n")
                  for (i in 1:dims$npred) {
                    cat(colnames(X)[i], "\n")
                  }

                  cat("\nResidual degrees of freedom = ", dims$rdf, "\n")
                  cat("Voxels = ", modData$nvox, "\n")
                  cat("FWHM = ", dims$fwhm, "\n")
                  cat("Resels = ", dims$resels, "\n\n")
                })
              )

#' @importFrom stats weights
#' @S3method weights rftModel
weights.rftModel <- function(object) {
  diag(object$W)
}
