#' Fit images to a linear model using random field theory
#'
#' @param modData
#' @param formula
#' @param weights vector of weights corresponding to each row of the design matrix and response matrix
#' @param checkMask logical. test to ensure only active voxels are tested
#' @param control a list of parameters for controlling the fitting process. (see \code{\link{rftControl}})
#' @param verbose activates verbose output
#'
#' @return An object of class "rftModel" see \code{\link{rftModel}} for more details
#'
#' @seealso \code{\link{modData}}
#' @description
#'
#' @References
#' Friston K.J., (1995) Statistical Parametric Maps in Functional Imaging: A General Linear Approach
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Stefan J.K., (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' @author Zachary P. Christensen
#' @examples
#'
#'
#' (fit <- mydata %>% rftLm(imat ~ var1) %>% summary(contrast = c(0, 1), cthresh = 100))
#'
#'
#' @export rftLm
rftLm <- function(imgData, formula, weights = NULL, control = list(), verbose = TRUE) {
  # get term specifying response/image group
  y <- all.vars(formula)[1]

  # isolate right hand side
  RHS <- delete.response(formula)
  X <- model.matrix(RHS, data = imgData$demog[imgData[y]$bool, ])
  cl <- match.call()


  fit <- rftFit(imgData, y, formula, X, K, V, weights, control, verbose)
  fit = fit$setResels(verbose = verbose)
  return(fit)
}

rftFit <- function(imgData, y, formula, X, K, V, weights, control,
                   verbose = NULL) {

  if (missing(K))
    K <- 1
  if (missing(xVi)) {
    xVi <- list()
    xVi$V <- diag(dataMod$nimg)
  }

  dims <- list(nimg = nrow(X), npred = ncol(X),
               nvox = ncol(imgData$imgList[y]$imageMatrix))
  fit <- rftModel(imgData = imgData, y = y, X = X, formula = formula, xVi = xVi, K = K, dims = dims, control = control, call = call)

  # get weight/whitening matrix: tcrossprod(W) = inv(V)
  if (missing(weights))
    fit = fit$setW()
  else
    fit = fit$setW(weights)

  # design space and projector matrix (pseudoinverse) for weighted least squares
  fit = fit$setKWX()
  fit = fit$setXX()

  # use non-sphericity xVi$V to compute effective degrees of freedom
  fit = fit$setXV()
  fit = fit$setDims()
  fit = fit$setBetaCov()

  fit = fit$setBeta()
  fit = fit$setResid()

  fit = fit$setCy(hsqr, UF, trMV)
  return(fit)
}

#'
#'
#'
#' @param method NA, ML, or REML are used to optimize fitting
#'
#' @export rftGlm
rftGlm <- function(imgData, formula, weights = NULL, method = "NA", control = list(),
                   verbose = TRUE) {
  # get term specifying response/image group
  y <- all.vars(formula)[1]

  # isolate right hand side
  RHS <- delete.response(formula)
  X <- model.matrix(RHS, data = imgData$demog[imgData[y]$bool, ])
  cl <- match.call()

  dims <- list(nimg = nrow(X), npred = ncol(X),
               nvox = ncol(imgData$imgList[y]$imageMatrix))

  fit <- rftModel(imgData = imgData, y = y, X = X, formula = formula,
     V = diag(dims$nimg), K = 1, dims = dims, control = control, call = call)

  # evoke ReML for hyperparameter estimation------------------------------------
  if (optim)
    fit <- rftModelOptimize(mod)

  # get weight/whitening matrix: tcrossprod(W) = inv(V)
  if (missing(W))
    fit = fit$setW()
  else
    fit = fit$setW(W)

  # design space and projector matrix (pseudoinverse) for weighted least squares
  fit = fit$setKWX()
  fit = fit$setXX()

  # use non-sphericity xVi$V to compute effective degrees of freedom
  fit = fit$setXV()
  fit = fit$setDims()
  fit = fit$setBetaCov()

  fit = fit$setBeta()
  fit = fit$setResid()
  fit = fit$setResels(verbose = verbose)
  return(fit)
}
