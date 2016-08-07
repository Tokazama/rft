# to do:
# model.matrix.iModel
# update iModel - still need to add updates methods for ...
# iModelRead - X and xVi

# iModel----
#' iModel class representing fitted models of image information
#' 
#' 
#'
#' coef fitted model.matrix residuals
#' 
#' @param X design matrix
#' @param y name of iGroup that corresponds to the response value
#' @param data iData object containing design information
#' @param weights optional vector or matrix of prior weights to be used in the fitting process
#' @param xVi 
#' @param control a list of parameters for controlling the fitting process. (see \code{link{iControl}} for more details)
#' @param filename
#' @param ...
#' 
#' 
#' 
#' @slot iData iData reference class object for fitted model
#' @slot y character naming the iGroup used within iData for the response to the fitted model
#' @slot X information about the design
#'  \item{X} {design matrix}
#'  \item{K} {filter information specific to each image}
#'  \item{W} {weights}
#'  \item{KWX} {weighted design matrix information}
#'  \item{XX} {pseudoinverse of X}
#'  \item{V} {non-sphericity matrix}
#'  \item{betaCov}
#'  \item{trRV} {trace of RV}
#'  \item{trRVRV} {trace of RVRV}
#'  \item{rdf} {residual degrees of freedom}
#' @slot beta beta coefficient matrix
#' @slot res residual matrix
#' @slot mrss mean residual sum of squares
#' @slot xVi information about intrinsic temporal non-sphericity
#'  \item{Vi} {}
#'  \item{V} {}
#'  \item{h} {hyperparameters}
#'  \item{Cy} {covariance of response matrix}
#' @slot dims model dimensions
#'  \item{npred} {number of predictors}
#'  \item{nimg} {number of images}
#'  \item{nvox} {number of voxels}
#'  \item{fwhm} {full-width at half-maxima}
#'  \item{resels} {resolution elements}
#'  \item{rpvImage} {resels per voxel image}
#' @slot call the original matched call information
#' @slot control control values for fitting the model (see \code{\link{iControl}})
#'
#' @export iModel
iModel <- setClass("iModel", 
                   slots = list(
                     iData = "iData",
                     y = "character",
                     X = "list",
                     beta = "DataSet",
                     res = "DataSet",
                     mrss = "DataSet",
                     xVi = "list",
                     dims = "list",
                     method = "character",
                     control = "list",
                     location = "character")
                   )

setMethod("initialize", "iModel", function(.Object, X, y = "unnamed", data, weights = NULL, xVi, control, filename, ...) {
  ll <- list(...)
  
  data <- if (missing(data)) iData()
  if (class(data) != "iData")
    stop("data must be an iData object")
  .Object@iData <- data
  if (length(y) > 1)
    stop("iModel currently only supports one response.")
  .Object@y <- y
  
  control <- iControl()
  if (!missing(control))
    control[names(control)] <- control
  .Object@control <- control
  
  # creae X list
  if(missing(X)) {
    .Object@X <- list()
    .Object@X$X <- matrix(1, 1, 1)
  }  else {
    .Object@X <- list()
    .Object@X$X <- X
  }
  
  # set dims
  dims <- list()
  dims$nimg <- nrow(.Object@X$X)
  dims$npred <- ncol(.Object@X$X)
  dims$nvox <- ncol(data@iList[[y]]@iMatrix[])
  if (control$rft) {
    dims$fwhm <- c()
    dims$resels <- c()
  }
  .Object@dims <- dims
  
  # set xVi
  if (missing(xVi)) {
    .Object@xVi <- list()
    .Object@xVi$V <- diag(dims$nimg)
  } else
    .Object@xVi <- xVi
  
  # set K
  K <- data@iList[[y]]@K
  if (class(K) == "numeric")
    .Object@K <- K
  else if (class(K) == "data.frame")
    .Object@K <- iFilter(K)
  
  # set X
  if (is.null(weights)) {
    iV <- sqrt(MASS::ginv(.Object@xVi$V))
    .Object@X$W <- iV * (abs(iV) > 1e-6)
  } else {
    if (class(weights) == "numeric" && length(weights) == dims$nimg)
      .Object@X$W <- diag(weights)
    else if (all(dim(weights) == dims$nimg))
      .Object@X$W <- weights
    else
      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
  }
  
  # set weighted design matrix
  .Object@X$KWX <- .setx(iFilter(.Object@X$K, .Object@X$W %*% .Object@X$X))
  # pseudoinverse of X
  .Object@X$XX <- .pinvx(.Object@X$KWX)
  .Object@X$V <- iFilter(.Object@X$K, iFilter(.Object@X$K, .Object@X$W %*% tcrossprod(.Object@xVi$V, .Object@X$W)))
  .Object@X$betaCov <- .Object@X$XX %*% tcrossprod(.Object@X$V, .Object@X$XX)
  out <- .trRV(.Object@X$KWX, .Object@X$V)
  .Object@X$trRV <- out$trRV
  .Object@X$trRVRV <- out$trRVRV
  .Object@X$rdf <- out$rdf
  
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  if (missing(filename)) {
    .Object@location <- tempfile(fileext = ".h5")
    file <- h5file(.Object@location)
  } else {
    if (file.exists(filename))
      stop("filename already exists.")
    .Object@location <- filename
    file <- h5file(filename)
  }
  
  
  if (is.function(.Object@control$chunksize))
    .Object@control$chunksize <- .Object@control$chunksize(.Object@dims$nimg)
  
  if (.Object@control$chunksize < 1)
    .Object@control$chunksize <- 1
  end <- .Object@control$chunksize
  nchunks <- floor(.Object@dims$nvox / end)
  start <- 1 + end
  if (.Object@dims$nvox > 1) {
    # setting chunksize for h5 data matrices
    if (is.null(ll$res)) {
      file["beta"] <- data.matrix(matrix(0, .Object@dims$nimg, end))
      file["beta"] <- cbind(file["beta"], data.matrix(matrix(0, .Object@dims$npred, .Object@dims$nvox - end)))
    } else
      file["beta"] <- ll$beta
    
    if (is.null(ll$res)) {
      file["res"] <- data.matrix(matrix(0, .Object@dims$nimg, end))
      file["res"] <- cbind(file["res"], data.matrix(matrix(0, .Object@dims$nimg, .Object@dims$nvox - end)))
      .Object@res <- file["res"]
    } else
      file["res"] <- ll$beta
    
    if (is.null(ll$mrss)) {
      file["mrss"] <- data.matrix(matrix(0, 1, end))
      file["mrss"] <- cbind(file["mrss"], data.matrix(matrix(0, 1, .Object@dims$nvox - end)))
    } else
      file["mrss"] <- ll$betax
  } else {
    file["beta"] <- matrix(0, .Object@dims$npred, .Object@dims$nvox)
    file["res"] <- matrix(0, .Object@dims$nimg, .Object@dims$nvox)
    file["mrss"] <- matrix(0, 1, .Object@dims$nvox)
  }
  .Object@mrss <- file["mrss"]
  .Object@beta <- file["beta"]
  .Object@res <- file["res"]
  
  return(.Object)
})

setMethod("show", "iModel", function(object) {
  cat("iModel object fit by call to ", object@method[1], " \n")
  cat("                  Predictors = ", object@dims$npred, "\n")
  cat("                    Names:\n")
  for (i in seq_len(object@dims$npred))
    cat("                   ", colnames(object@X$X)[i], "\n")
  cat("                      Images = ", object@dims$nimg, "\n")
  cat(" Residual degrees of freedom =", object@dims$rdf, "\n")
  cat("                      Voxels = ", object@dims$nvox, "\n")
  cat("                Optimization = ", object@control$opt, "\n")
})

#' @details get fitted coefficients from iModel
#' @describeIn iModel
setMethod("coef", "iModel", function(object) {
  object@beta[]
})

#' @details get fitted values from iModel
#' @describeIn iModel
setMethod("fitted", "iModel", function(object) {
  object@X$X %*% object@beta[]
})

#' @details get residuals from iModel object
#' @describeIn iModel
setMethod("resid", "iModel", function(object) {
  object@res[]
})

#' @details 
#' @describeIn iModel
iModelRead <- function(filename, iData_dirname) {
  if (!file.exists(filename))
    stop("filename does not exist.")
  x <- iModel()
  
  x@location <- filename
  x@y <- file["y"][]
  x@method <- file["method"][]
  
  # read control
  x@control$cf <- file["control/cf"][]
  x@control$mi <- file["control/mi"][]
  x@control$scr <- file["control/scr"][]
  x@control$sar <- file["control/sar"][] 
  x@control$tt <- file["control/tt"][]
  x@control$pval <- file["control/pval"] 
  x@control$pp <- file["control/pp"][] 
  x@control$n <- file["control/n"][]
  x@control$iso <- file["control/iso"][]  
  x@control$os <- file["control/os"][]
  x@control$rft <- file["control/rft"][] 
  x@control$opt <- file["control/opt"][] 
  x@control$chunksize <- file["control/chunksize"][] 
  
  # read dims
  x@dims$npred <- file["dims/npred"]
  x@dims$nvox <- file["dims/nvox"]
  x@dims$nimg <- file["dims/nimg"] 
  if (x@control$rft) {
    x$dims$fwhm <- file["dims/fwhm"]
    x@dims$resels <- file["dims/resels"]
    x@dims$rpvImage <- as.antsImage(file["dims/rpvImage"][])
    k = antsSetSpacing(x@dims$rpvImage, h5attr(file["dims/rpvImage"] , "spacing"))
    k = antsSetDirection(x@dims$rpvImage, h5attr(file["dims/rpvImage"] , "direction"))
    k = antsSetOrigin(x@dims$rpvImage, h5attr(file["dims/rpvImage"] , "origin"))
  }
  
  # read X
  x@X$X <- file["X/X"][]
  colnames(x@X$X) <- h5attr(file["X/X"], "colnames")[]
  x@X$W <- file["X/W"][]
  x@X$XX <- file["X/XX"][]
  x@X$V <- file["X/V"][]
  x@X$betaCov <- file["X/betaCov"][]
  x@X$trRV <- file["X/trRV"][]
  x@X$trRVRV <- file["X/trRVRV"][]
  x@X$rdf <- file["X/rdf"][]
  x@X$KWX$X <- file["X/KWX/X"][]
  x@X$KWX$v <- file["X/KWX/v"][]
  x@X$KWX$u <- file["X/KWX/u"][]
  x@X$KWX$d <- file["X/KWX/d"][]
  x@X$KWX$tol <- file["X/KWX/tol"][]
  x@X$KWX$rk <- file["X/KWX/rk"][]
  
  
  K <- data.frame(Filters = x@file[file.path("K", "Filters")][], 
                  HParam = x@file[file.path("K", "HParam")][],
                  RT = x@file[file.path("K", "RT")][])
  x@X$K <- iFilter(K)
  
  # read xVi
  if (file["xVi/Vi"][] != "null")
    x@xVi$Vi <- file["xVi/Vi"][]
  if (file["xVi/V"][] != "null")
    x@xVi$V <- file["xVi/V"][]
  if (file["xVi/h"][] != "null")
    x@xVi$h <- file["xVi/h"][]
  if (file["xVi/Cy"][] != "null")
    x@xVi$Cy <- file["xVi/Cy"][]
  
  return(x)
}

#' @details 
#' @describeIn iModel
iModelWrite <- function(x, filename, iData_dirname) {
  if (class(x) != "iModel")
    stop("x must be of class iModel.")
  if (!missing(iData_dirname))
    iDataWrite(x@iData, iData_dirname)
  if (file.exists(filename)) {
    if (x@location != filename)
      stop("filename already exists.")
  }
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  file.copy(x@location, filename)
  
  file["y"] <- x@y
  file["method"] <- x@method
  
  # write control
  file["control/cf"] <- x@control$cf
  file["control/mi"] <- x@control$mi
  file["control/scr"] <- x@control$scr
  file["control/sar"] <- x@control$sar
  file["control/tt"] <- x@control$tt
  file["control/pval"] <- x@control$pval
  file["control/pp"] <- x@control$pp
  file["control/n"] <- x@control$n
  file["control/iso"] <- x@control$iso 
  file["control/os"] <- x@control$os
  file["control/rft"] <- x@control$rft
  file["control/opt"] <- x@control$opt
  file["control/chunksize"] <- x@control$chunksize
  
  # write dims
  file["dims/npred"] <- x@dims$npred
  file["dims/nvox"] <- x@dims$nvox
  file["dims/nimg"] <- x@dims$nimg
  if (x@control$rft) {
    file["dims/fwhm"] <- x@dims$fwhm
    file["dims/resels"] <- x@dims$resels
    file["dims/rpvImage"] <- as.array(x@dims$rpvImage)
    h5attr(file["dims/rpvImage"] , "spacing") <- antsGetSpacing(x@dims$rpvImage)
    h5attr(file["dims/rpvImage"] , "direction") <- antsGetDirection(x@dims$rpvImage)
    h5attr(file["dims/rpvImage"] , "origin") <- antsGetOrigin(x@dims$rpvImage)
  }
  
  # write X
  file[file.path("K", "Filters")] <- x@iData@iList[[x@y]]@K$Filters
  file[file.path("K", "HParam")] <- x@iData@iList[[x@y]]@K$HParam
  file[file.path("K", "RT")] <- x@iData@iList[[x@y]]@K$RT
  
  file["X/X"] <- data.matrix(x@X$X)
  h5attr(file["X/X"], "colnames") <- colnames(x@X$X)
  file["X/W"] <- data.matrix(x@X$W)
  file["X/XX"] <- data.matrix(x@X$XX)
  file["X/V"] <- data.matrix(x@X$V)
  file["X/betaCov"] <- data.matrix(x@X$betaCov)
  file["X/trRV"] <- x@X$trRV
  file["X/trRVRV"] <- x@X$trRVRV
  file["X/rdf"] <- x@X$rdf
  file["X/KWX/X"] <- x@X$KWX$X
  file["X/KWX/v"] <- x@X$KWX$v
  file["X/KWX/u"] <- x@X$KWX$u
  file["X/KWX/d"] <- x@X$KWX$d
  file["X/KWX/tol"] <- x@X$KWX$tol
  file["X/KWX/rk"] <- x@X$KWX$rk
  
  # write xVi
  file["xVi/Vi"] <- if (is.null(x@xVi$Vi)) "null" else x@xVi$Vi
  file["xVi/V"] <- if (is.null(x@xVi$V)) "null" else x@xVi$V
  file["xVi/h"] <- if (is.null(x@xVi$h)) "null" else x@xVi$h
  file["xVi/Cy"] <- if (is.null(x@xVi$Cy)) "null" else x@xVi$Cy
  
  return(TRUE)
}

#' @details solve already initialized iModel for coefficients, residuals, and mean residual sum of squares
#' @describeIn iModel
iModelSolve <-  function(x, sr = TRUE, verbose = TRUE) {
  end <- x@control$chunksize
  nchunks <- floor(x@dims$nvox / end)
  start <- 1 + end
  for (j in seq_len(nchunks)) {
    if (j == nchunks)
      vrange <- start:x@dims$nvox
    else
      vrange <- start:end
    
    KWY <- iFilter(x@X$K, x@X$W %*% x@iData[[x@y]]@iMatrix[, vrange])
    x@beta[, vrange] <- x@X$XX %*% x@KWY
    x@res[, vrange] <- .res(X$KWX, KWY)
    x@mrss[, vrange] <- colSums(x@res[, vrange]^2) / x@X$trRV
    if (sr)
      x@res[, vrange] <- t(t(x@res[, vrange]) * (1 / as.numeric(x@mrss[, vrange])))
    
    start <- start + x@control$chunksize
    end <- end + x@control$chunksize
  }
  return(x)
}

#' @details primarily used to update X slot after optimization
#' @describeIn iModel
iModelUpdate <- function(x, ...) {
  ll <- c(...)
  if (!is.null(ll$weights)) {
    if (class(weights) == "numeric" && length(weights) == dims$nimg)
      x@X$W <- diag(weights)
    else if (all(dim(weights) == dims$nimg))
      x@X$W <- weights
    else
      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
  }
  
  # set weighted design matrix
  x@X$KWX <- .setx(iFilter(x@X$K, x@X$W %*% x@X$X))
  # pseudoinverse of X
  x@X$XX <- .pinvx(x@X$KWX)
  x@X$V <- iFilter(x@X$K, iFilter(x@X$K, x@X$W %*% tcrossprod(x@xVi$V, x@X$W)))
  x@X$betaCov <- x@X$XX %*% tcrossprod(x@X$V, x@X$XX)
  out <- .trRV(x@X$KWX, x@X$V)
  x@X$trRV <- out$trRV
  x@X$trRVRV <- out$trRVRV
  x@X$rdf <- out$rdf
  return(x)
}

#' Control parameters for RFT based analyses
#' 
#' @return 
#' \item{cf} {critical F-threshold for selecting voxels over which the non-sphericity is estimated (default = 0.001)}
#' \item{mi} {maximum iterations for optimizing fitted models}
#' \item{scr} {logical. scale residuals? (default = \code{TRUE})}
#' \item{sar} {number of residual images to sample for estimating the FWHM (default = \code{64})}
#' \item{tt} {threshType (see \code{\linke{statFieldThresh}})} 
#' \item{pval} {thresh p-value}
#' \item{pp} {primary/initial p-value threshold used for FDR thresholding}
#' \item{n} {images in conjunction (default = \code{1})}
#' \item{iso} {logical. should images be assumed to be isotropic? (default = \code{TRUE})}
#' \item{os} {offset weighting for iteratively reweighted least squares (default = \code{3})}
#' \item{rft} {logical. should voxels be estimated in resel space for random field theory analysis (default = \code{TRUE})}
#' \item{opt} {optimization method}
#' @export iControl
iControl <- function() {
  list(cf = 0.05,
       mi = 200,
       scr = TRUE,
       sar = 64,
       tt = "pRFT",
       pval = 0.05,
       pp = 0.01,
       n = 1,
       iso = TRUE,
       os = 3,
       rft = TRUE,
       opt = "none",
       chunksize = function(nimg) ((2^23) / nimg)
  )
}