#' Fit images to a generalized linear model using random field theory
#' 
#' @param images either an image matrix or image list
#' @param RH the right hand side of a formula (i.e. " ~ pred1 + pred2)
#' @param intercept should the design matrix include an intercept. logical (default = FALSE)
#' @param mask binary mask of class antsImage
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which rftGlm is called.
#' @param contrastMatrix a matrix in which each row represents a unique contrast and each column the weights for each predictor
#' @param fieldType the type of contrast to compute corresponding to each row of contrastMatrix ("F" for F-statistic & "T" for t-statistic)
#' @param statdir optional directory to save results
#' @param sample the amount of residual images used to estimate smoothness parameters (default uses all residuals)
#' @param verbose
#' @param smooth_sigma
#' @param filterfun filter/whiten function used to filter response and predictors prior to fitting (i.e. filterfun = whiten; default filterfun = NULL)
#' @param weights vector of weights corresponding to each row of the design matrix and response matrix
#' @param V optional non-sphericity matrix (if not specified this is estimated)
#' @param control a list of parameters for controlling the fitting process. (see rftControl)
#' 
#' @return
#' \item{call} {the matched call}
#' \item{dof} {degrees of freedom}
#' \item{weights} {weight matrix where the diagonal corresponds to each rows weight}
#' \item{coefficients} {model coefficients}
#' \item{rss} {residual sum of squares}
#' \item{residuals} {a matrix of residuals for each image or the directory where they were saved (if statdir was specified)}
#' \item{parameters}
#' \itemize{
#'   \item{S} {array of non-sphericity components}
#'   \item{V} {non-sphericity matrix}
#'   \item{h} {hyperparameters}
#'   \item{Cy} {covariance matrix of images}
#'   \item{R} {residual forming projection matrix}
#'   }
#' \item{volumes}
#' \itemize{
#'   \item{mask}
#'   \item{rpvImage}
#'   \item{fwhm}
#'   \item{resels}
#'   \item{voxels}
#'   \item{nimg}
#'   \item{npred}
#'   }
#' \item{contrasts}
#' \itemize{
#'   \item{contrastMatrix}
#'   \item{fieldType}
#' }
#' \item{statdir}
#' @details
#' 
#' 
#' Internal functions:
#' glFormula      rftModel
#' mkGlmerDevFun  rftREML
#' optimizeGlmer  rftREML
#' updateGlmer    update.rftModel
#' mkMerMod       rftModel
#' 
#' @author Zachary P. Christensen
#' 
#' 
#' 
#' 
#' 
#' 
#' @export rftGlm
rftGlm <- function(images, RH, intercept = FALSE, mask, data = NULL, contrastMatrix = NULL, statdir = NULL, sample = NULL, verbose = TRUE) {
  
  if (class(images) == "matrix") {
    n <- dim(images)[[1]]
    v <- dim(images)[[2]]
  } else {
    n <- length(images)
    v <- sum(mask)
  }
  
  
  # estimate parameters--------------------------------------------------------
  if (is.null(data))
    data <- environment(RH)
  X <- model.matrix(update(RH, ~ . - 1), data = data)
  
  # fit model------------------------------------------------------------------
  if (missing(weights)) {
    mymod <- rftModel(y, X, newmask, contrastMatrix, diag(n), findVar = TRUE, findResels = FALSE, "rftGlm")
    V <- rftREML(mymod$Cy, X, mymod$V)$V
    iV <- sqrt(MASS::ginv(V))
    weights <- iV * (abs(iV) > 1e-6)
  }
  mymod <- update(mymod, weights = weights, findResels = TRUE, findVar = FALSE)
  mymod$call <- call
  mymod
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
