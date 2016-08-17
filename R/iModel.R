# TO DO:
# replace initialize with iModelMake
# plot.iModel
#
# Error in iModelMake(X = z$X, y = z$y[i], iData = z$iData) : 
# trying to get slot "K" from an object of a basic class ("NULL") with no slots
#

#' iModel class representing fitted models of image information
#' 
#' 
#'
#' 
#' @param X design matrix
#' @param y name of iGroup that corresponds to the response value
#' @param data iData object containing design information
#' @param weights optional vector or matrix of prior weights to be used in the fitting process
#' @param xVi sphericity measurements
#' @param control a list of parameters for controlling the fitting process. (see \code{\link{iControl}} for more details)
#' @param filename optional filename location to save object to 
#' @param ... named arguments
#' 
#' @slot iData iData reference class object for fitted model
#' @slot y character naming the iGroup used within iData for the response to the fitted model
#' @slot X information about the design
#' \itemize{
#'  \item{X} {design matrix}
#'  \item{K} {filter information specific to each image}
#'  \item{W} {weights}
#'  \item{KWX} {weighted design matrix information}
#'  \item{pKX} {pseudoinverse of X}
#'  \item{V} {non-sphericity matrix}
#'  \item{betaCov}
#'  \item{trRV} {trace of RV}
#'  \item{trRVRV} {trace of RVRV}
#'  \item{rdf} {residual degrees of freedom}
#' }
#' @slot beta beta coefficient matrix
#' @slot res residual matrix
#' @slot mrss mean residual sum of squares
#' @slot xVi information about intrinsic temporal non-sphericity
#' \itemize{
#'  \item{Vi} {list of non-sphericity components}
#'  \item{V} {non-sphericity matrix (Cov(e) = sigma^2*V)}
#'  \item{h} {hyperparameters}
#'  \item{Cy} {covariance of response matrix}
#' }
#' @slot dims model dimensions
#' \itemize{
#'  \item{npred} {number of predictors}
#'  \item{nimg} {number of images}
#'  \item{nvox} {number of voxels}
#'  \item{fwhm} {full-width at half-maxima}
#'  \item{resels} {resolution elements}
#'  \item{rpvImage} {resels per voxel image}
#' }
#' @slot call the original matched call information
#' @slot C list of contrasts
#' \itemize{
#'  \item{name} {name of the contrast}
#'  \item{c} {contrast weights}
#'  \item{X1} {Reamaining design space (orthogonal to X0)}
#'  \item{X0} {Reduced design matrix.}
#'  \item{iX0} {Indicates how contrat was specified}
#'  \item{dims} {dimensions for contrast}
#'  \itemize{
#'    \item{trRV} {trace of RV}
#'    \item{trRVRV} {trace of RVRV}
#'    \item{idf} {degrees of interest}
#'  }
#'  \item{fieldType} {type of statyistical field being fitted}
#'  \item{contrastImage} {object of class antsImage representing contrast}
#'  \item{clusterImage} {object of class antsImage representing the thresholded contrastImage}
#'  \item{results} {either a list or data.frame containing the results of a given contrast}
#'  \item{sthresh} {statistical threshold}
#'  \item{cthresh} {cluster threshold}
#' }
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
                     C = "list",
                     method = "character",
                     control = "list",
                     location = "character")
                   )

iModelMake <- function(X, y, iData, weights = NULL, xVi, control, filename, ...) {
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  out <- iModel()
  
  ll <- list(...)
  
  if(missing(iData))
    out@iData <- iData()
  else 
    out@iData <- iData
  
  if (length(y) > 1)
    stop("iModel currently only supports one response.")
  if (missing(y))
    out@y <- out@iData@iList[[1]]@name
  else
    out@y <- y

  control <- iControl()
  if (!missing(control))
    control[names(control)] <- control
  out@control <- control
  
  # create X list
  if(missing(X)) {
    out@X <- list()
    out@X$X <- matrix(1, 1, 1)
  }  else {
    out@X <- list()
    out@X$X <- X
  }
  
  # set dims
  dims <- list()
  dims$nimg <- nrow(out@X$X)
  dims$npred <- ncol(out@X$X)
  dims$nvox <- ncol(out@iData@iList[[out@y]])
  if (control$rft) {
    dims$fwhm <- c()
    dims$resels <- c()
  }
  out@dims <- dims
  
  # set xVi
  if (missing(xVi))
    out@xVi$V <- diag(dims$nimg)
  else
    out@xVi <- xVi
  
  # set K
  K <- out@iData@iList[[out@y]]@K
  if (class(K) == "numeric")
    out@X$K <- K
  else if (class(K) == "data.frame")
    out@X$K <- iFilter(K)
  
  # set X
  if (is.null(weights)) {
    iV <- sqrt(MASS::ginv(out@xVi$V))
    out@X$W <- iV * (abs(iV) > 1e-6)
  } else {
    if (class(weights) == "numeric" && length(weights) == dims$nimg)
      out@X$W <- diag(weights)
    else if (all(dim(weights) == dims$nimg))
      out@X$W <- weights
    else
      stop("weights must be a matrix of nimg x nimg or a vector of length nimg")
  }
  
  # set weighted design matrix
  out@X$KWX <- .setx(iFilter(out@X$K, out@X$W %*% out@X$X))
  # pseudoinverse of X
  out@X$pKX <- .pinvx(out@X$KWX)
  out@X$V <- iFilter(out@X$K, iFilter(out@X$K, out@X$W %*% tcrossprod(out@xVi$V, out@X$W)))
  out@X$betaCov <- out@X$pKX %*% tcrossprod(out@X$V, out@X$pKX)
  out@X[c("trRV", "trRVRV", "rdf")] <- .trRV(out@X$KWX, out@X$V)
  
  if (!usePkg("h5"))
    stop("Please install package h5 in order to use this function.")
  if (missing(filename)) {
    out@location <- tempfile(fileext = ".h5")
    file <- h5file(out@location)
  } else {
    if (file.exists(filename))
      stop("filename already exists.")
    out@location <- filename
    file <- h5file(filename)
  }
  
  end <- out@iData@iList[[out@y]]@iMatrix@chunksize[2]
  nchunks <- floor(out@dims$nvox / end)
  if (nchunks > 1) {
    # setting chunksize for h5 data matrices
    if (is.null(ll$res)) {
      file["beta"] <- data.matrix(matrix(0, out@dims$npred, end))
      file["beta"] <- cbind(file["beta"], data.matrix(matrix(0, out@dims$npred, out@dims$nvox - end)))
    } else
      file["beta"] <- ll$beta
    
    if (is.null(ll$res)) {
      file["res"] <- data.matrix(matrix(0, out@dims$nimg, end))
      file["res"] <- cbind(file["res"], data.matrix(matrix(0, out@dims$nimg, out@dims$nvox - end)))
      out@res <- file["res"]
    } else
      file["res"] <- ll$beta
    
    if (is.null(ll$mrss)) {
      file["mrss"] <- data.matrix(matrix(0, 1, end))
      file["mrss"] <- cbind(file["mrss"], data.matrix(matrix(0, 1, out@dims$nvox - end)))
    } else
      file["mrss"] <- ll$betax
  } else {
    file["beta"] <- matrix(0, out@dims$npred, out@dims$nvox)
    file["res"] <- matrix(0, out@dims$nimg, out@dims$nvox)
    file["mrss"] <- matrix(0, 1, out@dims$nvox)
  }
  out@mrss <- file["mrss"]
  out@beta <- file["beta"]
  out@res <- file["res"]
  
  return(out)
}

#' @export
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
  
  if (object@control$rft && length(object@dims$fwhm > 0)) {
    cat("                        FWHM = ", round(object@dims$fwhm, 2), "\n")
    cat("                      Resels = ", round(object@dims$resels), "\n")
  }
  
  ncon <- length(object@C)
  if (ncon == 0)
    cat("Contrasts not set.\n")
  else {
    for (i in seq_len(ncon)) {
      cat("Contrast", x@C[[i]]$name, "\n")
      cat(" Contrast weights: \n")
      print(t(object@c))
      
      if (object@control$rft) {
        cat(" Set-level: \n")
        cat(" Clusters = ", ncol(object@C[[i]]$results$clusterLevel), "\n")
        cat(" p-value = ", object@C[[i]]$results$setLevel, "\n\n")
        
        cat("  Cluster-Level: \n")
        print(round(object@C[[i]]$results$clusterLevel, 3))
        cat("\n")
        
        cat(" Peak-level: \n")
        print(round(object@C[[i]]$results$peakLevel, 3))
      } else {
        print(object@C[[i]]$results)
      }
      cat("Interest degrees of freedom = ", object@C[[i]]$dims$idf, "\n")
      cat("Statistical threshold = ", round(object@C[[i]]$sthresh, 2), "\n")
      cat("Cluster threshold = ", object@C[[i]]$cthresh, "\n")
      cat("\n___\n\n")
    }
  }
})

#' iModel Methods
#' 
#' @param x,object Object of class iModel.
#' @param filename h5 file to save iModel to.
#' @param iData_dirname Directory for iData component of iModel.
#' @param verbose enables verbose output. (default = \code{TRUE})
#' @param ... additional named arguments passed to \code{iModelUpdate}
#' @param contrastMatrix matrix of contrasts
#' @param cthresh cluster level threshold
#' @param fieldType statistical field type for contrast
#' @name iModel-methods
NULL

#' @export
#' @docType methods
#' @details \strong{coef} Retrieve fitted coefficients from iModel object.
#' @rdname iModel-methods
setMethod("coef", "iModel", function(object) {
  object@beta[]
})

#' @export
#' @docType methods
#' @details \strong{fitted} Retrieve fitted values from iModel object.
#' @rdname iModel-methods
setMethod("fitted", "iModel", function(object) {
  object@X$X %*% object@beta[]
})

#' @export
#' @docType methods
#' @details \strong{resid} Retrieve residuals from iModel object.
#' @rdname iModel-methods
setMethod("resid", "iModel", function(object) {
  object@res[]
})

#' @export
#' @docType methods
#' @details \strong{iModelRead} read/load iModel object
#' @rdname iModel-methods
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
  x@X$pKX <- file["X/pKX"][]
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
  
  # read contrasts
  ncon <- file["C/count"][]
  for (i in seq_len(x@C)) {
    cname <- paste("C/C", i, sep = "")
    x@C[[i]]$name <- file[file.path(cname, "name")][]
    x@C[[i]]$c <- file[file.path(cname, "c")][]
    x@C[[i]]$X1 <- file[file.path(cname, "X1")][] 
    x@C[[i]]$X0 <- file[file.path(cname, "X0")][] 
    x@C[[i]]$iX0file[file.path(cname, "iX0")][] 
    x@C[[i]]$fieldType <- file[file.path(cname, "fieldType")][]
    x@C[[i]]$results <- file[file.path(cname, "results")][]
    x@C[[i]]$sthresh <- file[file.path(cname, "sthresh")][]
    x@C[[i]]$cthresh <- file[file.path(cname, "cthresh")][]
    
    x@C[[i]]$contrastImage <- as.antsImage(file[file.path(cname, "contrastImage")][])
    x@C[[i]]$contrastImage <- antsSetSpacingh5attr(file[file.path(cname, "contrastImage")], "spacing")
    x@C[[i]]$contrastImage <- antsSetDirection(h5attr(file[file.path(cname, "contrastImage")], "direction"))
    x@C[[i]]$contrastImage <- antsSetOrigin(h5attr(file[file.path(cname, "contrastImage")], "origin"))
    
    x@C[[i]]$clusterImage <- as.antsImage(file[file.path(cname, "clusterImage")][])
    x@C[[i]]$clusterImage <- antsSetSpacingh5attr(file[file.path(cname, "clusterImage")], "spacing")
    x@C[[i]]$clusterImage <- antsSetDirection(h5attr(file[file.path(cname, "clusterImage")], "direction"))
    x@C[[i]]$clusterImage <- antsSetOrigin(h5attr(file[file.path(cname, "clusterImage")], "origin"))
    
    x@C[[i]]$dims$trMV <- file[file.path(cname, "dims", "trMV")][]
    x@C[[i]]$dims$trMVMV <- file[file.path(cname, "dims", "trMVMV")][]
    x@C[[i]]$dims$idf <- file[file.path(cname, "dims", "idf")][]
  }
  return(x)
}

#' @export
#' @docType methods
#' @details \strong{iModelWrite} read/load iModel object
#' @rdname iModel-methods
# @describeIn iModel write/save iModel object
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
  file["X/pKX"] <- data.matrix(x@X$pKX)
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
  
  # write contrasts
  file["C/count"] <- length(x@C)
  for (i in seq_len(x@C)) {
    cname <- paste("C/C", i, sep = "")
    file[file.path(cname, "name")] <- x@C[[i]]$name
    file[file.path(cname, "c")] <- x@C[[i]]$c
    file[file.path(cname, "X1")] <- x@C[[i]]$X1
    file[file.path(cname, "X0")] <- x@C[[i]]$X0
    file[file.path(cname, "iX0")] <- x@C[[i]]$iX0
    file[file.path(cname, "fieldType")] <- x@C[[i]]$fieldType
    file[file.path(cname, "results")] <- x@C[[i]]$results
    file[file.path(cname, "sthresh")] <- x@C[[i]]$sthresh
    file[file.path(cname, "cthresh")] <- x@C[[i]]$cthresh
    
    file[file.path(cname, "contrastImage")] <- as.array(x@C[[i]]$contrastImage)
    h5attr(file[file.path(cname, "contrastImage")], "spacing") <- antsGetSpacing(x@C[[i]]$contrastImage)
    h5attr(file[file.path(cname, "contrastImage")], "direction") <- antsGetDirection(x@C[[i]]$contrastImage)
    h5attr(file[file.path(cname, "contrastImage")], "origin") <- antsGetOrigin(x@C[[i]]$contrastImage)
    
    file[file.path(cname, "clusterImage")] <- as.array(x@C[[i]]$clusterImage)
    h5attr(file[file.path(cname, "clusterImage")], "spacing") <- antsGetSpacing(x@C[[i]]$clusterImage)
    h5attr(file[file.path(cname, "clusterImage")], "direction") <- antsGetDirection(x@C[[i]]$clusterImage)
    h5attr(file[file.path(cname, "clusterImage")], "origin") <- antsGetOrigin(x@C[[i]]$clusterImage)
    
    file[file.path(cname, "dims", "trMV")] <- x@C[[i]]$dims$trMV
    file[file.path(cname, "dims", "trMVMV")] <- x@C[[i]]$dims$trMVMV
    file[file.path(cname, "dims", "idf")] <- x@C[[i]]$dims$idf
  }
  return(TRUE)
}

#' @export
#' @docType methods
#' @details \strong{iModelSolve} Solve already initialized iModel for coefficients, residuals, and mean residual sum of squares
#' @rdname iModel-methods
iModelSolve <-  function(x, verbose = TRUE) {
  end <- x@iData@iList[[x@y]]@iMatrix@chunksize
  nchunks <- floor(x@dims$nvox / end)
  start <- 1 + end
  for (j in seq_len(nchunks)) {
    if (j == nchunks)
      vrange <- start:x@dims$nvox
    else
      vrange <- start:end
    
    KWY <- iFilter(x@X$K, x@X$W %*% x@iData[[x@y]]@iMatrix[, vrange])
    x@beta[, vrange] <- x@X$pKX %*% x@KWY
    x@res[, vrange] <- .res(X$KWX, KWY)
    x@mrss[, vrange] <- colSums(x@res[, vrange]^2) / x@X$trRV
    if (x@control$scr)
      x@res[, vrange] <- t(t(x@res[, vrange]) * (1 / as.numeric(x@mrss[, vrange])))
    
    start <- start + x@control$chunksize
    end <- end + x@control$chunksize
  }
  return(x)
}

#' @export
#' @docType methods
#' @details \strong{iModelUpdate} Primarily used to update X slot after optimization
#' @rdname iModel-methods
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
  x@X$pKX <- .pinvx(x@X$KWX)
  x@X$V <- iFilter(x@X$K, iFilter(x@X$K, x@X$W %*% tcrossprod(x@xVi$V, x@X$W)))
  x@X$betaCov <- x@X$pKX %*% tcrossprod(x@X$V, x@X$pKX)
  out <- .trRV(x@X$KWX, x@X$V)
  x@X$trRV <- out$trRV
  x@X$trRVRV <- out$trRVRV
  x@X$rdf <- out$rdf
  return(x)
}

#' @export
#' @docType methods
#' @details \strong{model.matrix} Retrieve design matrix from iModel object.
#' @rdname iModel-methods
setMethod("model.matrix", "iModel", function(object) {
  return(object@X$X)
})

#' @export
#' @docType methods
#' @description \strong{summary} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-methods
setMethod("summary", "iModel", function(object, contrastMatrix, cthresh = 150, fieldType) {
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  
  if (length(contrast) != object@dims$npred)
    stop("the contrast length must be equal to the number of predictors")
  mask <- object@iData@iList[[object@y]]@mask
  
  # ensure that if previous contrasts were set they won't be erased
  start <- length(object@C) + 1
  ncon <- nrow(contrastMatrix) + start
  
  for (i in start:ncon) {
    contrast <- matrix(contrastMatrix[i, ], ncol = 1)
    tmp <- matrix(con, ncol = 1)
    rownames(tmp) <- colnames(object@X$X)
    c <- tmp
    object@C[[i]] <- .setcon(connames[i], fieldType[i], action, contrast, object@X$KWX)
    
    # X1 <<- .X1(.self, iModel$X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)
    names(object@C[[i]]$dims) <- c("trMV", "trMVMV", "idf")
    
    # solve contrast
    if (object@C[[i]]$fieldType == "T") {
      Vc <- crossprod(object@C[[i]]$c, object@covBeta) %*% object@C[[i]]$c
      se <- sqrt(object@mrss * as.numeric(Vc))
      tvec <- crossprod(object@C[[i]]$c, object@beta) / se
      object@C[[i]]$contrastImage <- makeImage(mask, tvec)
    } else if (object@C[[i]]$fieldType == "F") {
      h <- .hsqr(object@C[[i]], object@X$KWX)
      ss <- (rowSums((h %*% object@beta)^2) / object@C[[i]]$dims$trMV)
      object@C[[i]]$contrastImage <- makeImage(mask, ss / object@mrss)
    }
    
    if (!object@control$rft) {
      if (!object@control$iso)
        isotropic <- object$rpvImage
      else
        isotropoic <- NULL
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm,
                      c(object@C[[i]]$dims$idf, object$dims$edf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = cthresh[i], ll$threshType, ll$pval, ll$pp, ll$n, ll$statdir, ll$verbose)
      
      object@C[[i]]$sthresh <- z$threshold
      object@C[[i]]$clusterImage <- z$clusterImage
      object@C[[i]]$results$setLevel <- z$setLevel
      object@C[[i]]$results$peakLevel <- z$peakLevel
      object@C[[i]]$results$clusterLevel <- z$clusterLevel
    } else {
      # don't use random field theory for summary
      object@C[[i]]$sthresh <-
        statFieldThresh(object@C[[i]]$contrastImage, ll$pval, object@dims$nvox,
                        ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)),
                        c(object@C[[i]]$dims$idf, object$dims$rdf),
                        object@C[[i]]$fieldType, ll$tt, ll$pp)
      
      object@C[[i]]$clusterImage <-
        labelClusters(object@C[[i]]$contrastImage, minClusterSize = cthresh[i],
                      minThresh = object@C[[i]]$sthresh, maxThresh = Inf)
      
      object@C[[i]]$results <- labelStats(object@C[[i]]$contrastImage, object@C[[i]]$clusterImage)
    }
  }
  return(object)
})

#' Control parameters for RFT based analyses
#' 
#' @param cf critical F-threshold for selecting voxels over which the non-sphericity is estimated (default = \code{0.001})
#' @param mi maximum iterations for optimizing fitted models
#' @param scr logical. scale residuals? (default = \code{TRUE})
#' @param sar number of residual images to sample for estimating the FWHM (default = \code{64})
#' @param tt threshType (see \code{statFieldThresh})
#' @param pval thresh p-value
#' @param pp primary/initial p-value threshold used for FDR thresholding
#' @param n images in conjunction (default = \code{1})
#' @param iso logical. should images be assumed to be isotropic? (default = \code{TRUE})
#' @param os offset weighting for iteratively reweighted least squares (default = \code{3})
#' @param rft logical. should voxels be estimated in resel space for random field theory analysis (default = \code{TRUE})
#' @return 
#' 
#' @export iControl
iControl <- function(cf = 0.05, mi = 200, scr = TRUE, sar = 64, tt = "pRFT",
                     pval = 0.05, pp = 0.01, n = 1, iso = TRUE, os = 3, rft = TRUE) {
  list(cf = cf,
       mi =  mi,
       scr = scr,
       sar = sar,
       tt = tt,
       pval = pval,
       pp = pp,
       n = n,
       iso = iso,
       os = os,
       rft = rft)
}