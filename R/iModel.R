# TO DO:
# plot.iModel
#

#' Class iModel
#' 
#' Object containing a fitted statistical model for image data.
#' 
#' @slot iData iData reference class object for fitted model.
#' @slot y Character naming the iGroup used within iData for the response to the fitted model.
#' @slot X Information about the design:
#' \itemize{
#'  \item{X} {Design matrix.}
#'  \item{K} {Filter information specific to each image.}
#'  \item{W} {weights}
#'  \item{KWX} {Weighted design matrix information.}
#'  \item{pKX} {Pseudoinverse of X.}
#'  \item{V} {Non-sphericity matrix.}
#'  \item{betaCov} {Covariance of correlation coefficients.}
#'  \item{trRV} {Trace of RV.}
#'  \item{trRVRV} {Trace of RVRV.}
#'  \item{rdf} {Residual degrees of freedom.}
#' }
#' @slot beta Pointer beta matrix stored as h5 dataset.
#' @slot res Pointer residual matrix stored as h5 dataset.
#' @slot mrss Pointer mean residual sum of squares matrix stored as h5 dataset.
#' @slot xVi Information about intrinsic temporal non-sphericity:
#' \itemize{
#'  \item{Vi} {List of non-sphericity components.}
#'  \item{V} {Non-sphericity matrix (Cov(e) = sigma^2*V)}
#'  \item{h} {Hyperparameters.}
#'  \item{Cy} {Covariance of response matrix.}
#' }
#' @slot dims Model dimensions:
#' \itemize{
#'  \item{npred} {Number of predictors.}
#'  \item{nimg} {Number of images.}
#'  \item{nvox} {Number of voxels.}
#'  \item{fwhm} {full-width at half-maxima}
#'  \item{resels} {Resolution elements.}
#'  \item{rpvImage} {Resels per voxel image.}
#' }
#' @slot C list of contrasts:
#' \itemize{
#'  \item{name} {Name of the contrast.}
#'  \item{c} {Contrast weights.}
#'  \item{X1} {Reamaining design space (orthogonal to X0)}
#'  \item{X0} {Reduced design matrix.}
#'  \item{iX0} {Indicates how contrat was specified.}
#'  \item{dims} {Dimensions for contrast}
#'  \itemize{
#'    \item{trRV} {Trace of RV.}
#'    \item{trRVRV} {Trace of RVRV.}
#'    \item{idf} {Degrees of interest.}
#'  }
#'  \item{fieldType} {Type of statyistical field being fitted.}
#'  \item{contrastImage} {Object of class antsImage representing contrast.}
#'  \item{clusterImage} {Object of class antsImage representing the thresholded contrastImage.}
#'  \item{results} {Either a list or data.frame containing the results of a given contrast.}
#'  \item{sthresh} {Statistical threshold.}
#'  \item{cthresh} {Cluster threshold.}
#' }
#' @slot control Control values for fitting the model (see \code{\link{iControl}}).
#' 
#' @author Zachary P. Christensen
#' 
#' @seealso \code{\link{ilm}}
#'
#' @rdname iModel-class
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

#' @export 
#' @docType methods
#' @details \strong{iModelMake} Make iModel object.
#' @rdname iModel-method
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
    out@X$K <- .filter(K)
  
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
  out@X$KWX <- .setx(.filter(out@X$K, out@X$W %*% out@X$X))
  # pseudoinverse of X
  out@X$pKX <- .pinvx(out@X$KWX)
  out@X$V <- .filter(out@X$K, .filter(out@X$K, out@X$W %*% tcrossprod(out@xVi$V, out@X$W)))
  out@X$betaCov <- out@X$pKX %*% tcrossprod(out@X$V, out@X$pKX)
  tmp <- .trRV(out@X$KWX, out@X$V)
  out@X$trRV <- tmp$trRV
  out@X$trRVRV <- tmp$trRVRV
  out@X$rdf <- tmp$rdf
  
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
  
  cat("                  Predictors =", colnames(object@X$X), "\n")
  cat("                      Images = ", object@dims$nimg, "\n")
  cat(" Residual degrees of freedom = ", object@X$rdf, "\n")
  cat("                      Voxels = ", object@dims$nvox, "\n")
  cat("                Optimization = ", object@control$opt, "\n")
  cat("                    location = ", object@location, "\n")
  
  if (object@control$rft && length(object@dims$fwhm > 0)) {
    cat("                        FWHM = ", round(object@dims$fwhm, 2), "\n")
    cat("                      Resels = ", round(object@dims$resels), "\n")
  }
  
  cat("--- \n\n")
  
  ncon <- length(object@C)
  if (ncon == 0)
    cat("Contrasts not set.\n")
  else {
    for (i in seq_len(ncon)) {
      cat("Contrast", object@C[[i]]$name, "\n")
      cat(" Contrast weights: \n")
      c <- t(object@C[[i]]$c)
      colnames(c) <- colnames(object@X$X)
      print(c)
      cat("\n")
      if (object@control$rft) {
        if (class(object@C[[i]]$results) == "character") {
          cat(object@C[[i]]$results, "\n\n")
        } else {
          cat(" Set-level: \n")
          cat("  Clusters = ", ncol(object@C[[i]]$results$clusterLevel), "\n")
          cat("  p-value  = ", object@C[[i]]$results$setLevel, "\n\n")
          
          cat(" Cluster-Level: \n")
          print(round(object@C[[i]]$results$clusterLevel, 3))
          cat("\n")
          
          cat(" Peak-level: \n")
          print(round(object@C[[i]]$results$peakLevel, 3))
          cat("\n")
        }
      } else {
        if (class(object@C[[i]]$results) == "character") {
          cat(object@C[[i]]$results, "\n\n")
        } else
          print(object@C[[i]]$results)
      }
      cat("Interest degrees of freedom = ", object@C[[i]]$dims$idf, "\n")
      cat("      Statistical threshold = ", round(object@C[[i]]$sthresh, 2), "\n")
      cat("          Cluster threshold = ", object@C[[i]]$cthresh, "\n")
      cat("             Threshold type = ", object@C[[i]]$threshType, "\n")
      cat("---\n\n")
    }
  }
})

#' iModel Methods
#' 
#' @param x,object Object of class iModel.
#' @param filename h5 file to save iModel to.
#' @param iData_dirname Directory for iData component of iModel.
#' @param docname Prefix name for report documents.
#' @param imgdir Optional directory to save images to.
#' @param output_format File format for report to be rendered in (see rmarkdown::render).
#' @param verbose Enables verbose output. (default = \code{TRUE}).
#' @param ... Additional named arguments passed to \code{iModelUpdate}.
#' 
#' @author Zachary P. Christensen
#' 
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
#' @details \strong{iModelRead} Read/load iModel object.
#' @rdname iModel-methods
iModelRead <- function(filename, iData_dirname) {
  if (!file.exists(filename))
    stop("filename does not exist.")
  x <- iModel()
  file <- h5file(filename)
  
  # load big matrices
  x@mrss <- file["mrss"]
  x@beta <- file["beta"]
  x@res <- file["res"]
  
  x@location <- filename
  x@y <- file["y"][]
  x@method <- file["method"][]
  
  # read control
  x@control$cf <- file["control/cf"][]
  x@control$scr <- file["control/scr"][]
  x@control$sar <- file["control/sar"][] 
  x@control$pval <- file["control/pval"] 
  x@control$n <- file["control/n"][]
  x@control$iso <- file["control/iso"][]  
  x@control$os <- file["control/os"][]
  x@control$rft <- file["control/rft"][] 

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
  x@X$K <- .filter(K)
  
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
#' @details \strong{iModelWrite} Read/load iModel object.
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
  
  file <- h5file(filename)
  
  #### write big matrices
  chunksize <- x@beta@chunksize[2]
  chunkseq <- seq_len(chunksize)
  nvox <- x@iData@iList[[x@y]]@dim[2]
  nchunk <- floor(nvox / chunksize)
  
  # initialize files
  file["beta"] <- x@beta[, chunkseq]
  file["res"]  <- x@res[, chunkseq]
  file["mrss"] <- x@mrss[, chunkseq]
  
  beta <- file["beta"]
  res  <- file["res"]
  mrss <- file["mrss"]
  
  # starting on 2nd chunk fill in matrices
  for (i in seq_len(nchunk - 1)) {
    chunkseq <- chunkseq + chunksize
    beta <- cbind(beta, x@beta[, chunkseq])
    res  <- cbind(res, x@res[, chunkseq])
    mrss <- cbind(mrss, x@mrss[, chunkseq])
  }
  
  # fill in any remaining chunkage
  if (((nvox / chunksize) - nchunk) > 0) {
    beta <- cbind(beta, x@beta[, chunkseq[chunksize]:nvox])
    res  <- cbind(res, x@res[, chunkseq[chunksize]:nvox])
    mrss <- cbind(mrss, x@mrss[, chunkseq[chunksize]:nvox])
  }
  #### end big matrices
  
  
  file["y"] <- x@y
  file["method"] <- x@method
  
  # write control
  file["control/cf"] <- x@control$cf
  file["control/scr"] <- x@control$scr
  file["control/sar"] <- x@control$sar
  file["control/n"] <- x@control$n
  file["control/iso"] <- x@control$iso 
  file["control/os"] <- x@control$os
  file["control/rft"] <- x@control$rft

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
#' @details \strong{iModelSolve} Solve already initialized iModel for
#'  coefficients, residuals, and mean residual sum of squares.
#' @rdname iModel-methods
iModelSolve <-  function(x, verbose = TRUE) {
  chunksize <- x@iData@iList[[x@y]]@iMatrix@chunksize[2]
  nchunks <- floor(x@dims$nvox / chunksize)
  vrange <- seq_len(chunksize)
  if (verbose)
    progress <- txtProgressBar(min = 0, max = nchunks, style = 3)
  for (j in seq_len(nchunks)) {
    if (vrange[chunksize] > x@dims$nvox)
      vrange <- vrange[1]:x@dims$nvox
    
    KWY <- .filter(x@X$K, x@X$W %*% x@iData@iList[[x@y]]@iMatrix[, vrange])
    x@beta[, vrange] <- x@X$pKX %*% KWY
    x@res[, vrange] <- .res(x@X$KWX, KWY)
    x@mrss[, vrange] <- colSums(x@res[, vrange]^2) / x@X$trRV
    
    if (x@control$scr)
      x@res[, vrange] <- t(t(x@res[, vrange]) * (1 / as.numeric(x@mrss[, vrange])))
    
    vrange <- vrange + chunksize
    if (verbose)
      setTxtProgressBar(progress, j)
  }
  if (verbose)
    close(progress)
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
  x@X$KWX <- .setx(.filter(x@X$K, x@X$W %*% x@X$X))
  # pseudoinverse of X
  x@X$pKX <- .pinvx(x@X$KWX)
  x@X$V <- .filter(x@X$K, .filter(x@X$K, x@X$W %*% tcrossprod(x@xVi$V, x@X$W)))
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
#' @details \strong{report} Creates a .pdf and .html report of for iModel
#'  objects (requires package rmarkdown).
#' @rdname iModel-methods
report <- function(x, docname, surfimg) {
  if (!usePkg("rmarkdown"))
    stop("Please install package rmarkdown in order to use this function.")
  ncon <- length(x@C)
  if (missing(surfimg))
    surfimg <- x@iData@iList[[x@y]]@mask
  ## render images----
  for (i in seq_len(ncon)) {
    if (class(x@C[[i]]$results) != "character") {
      brain <- renderSurfaceFunction(surfimg = list(surfimg), funcimg = list(x@C[[i]]$clusterImage), alphasurf = 0.1, smoothsval = 1.5)
      tmp <- make3ViewPNG(fnprefix = paste(docname, "_", x@C[[i]]$name, sep = ""))
    }
  }
  
  md <- paste(docname, ".Rmd", sep =)
  zz <- file(md, open = "wt")
  sink(zz)
  sink(zz, type = "message")
  
  # describe the model----
  cat("## iModel object fit by call to ", x@method[1])
  if (x@control$rft)
    cat(" using random field theory. \n")
  else
    cat(". \n\n")
  
  cat("                Predictors = ", colnames(x@X$X), "\n\n")
  cat("                    Images = ", x@dims$nimg, "\n\n")
  cat("        Degrees of freedom = ", x@X$rdf, "\n\n")
  cat("                    Voxels = ", x@dims$nvox, "\n\n")
  cat("              Optimization = ", x@control$opt, "\n\n")
  
  if (x@control$rft && length(x@dims$fwhm > 0)) {
    cat("                      FWHM = ", round(x@dims$fwhm, 2), "\n\n")
    cat("                    Resels = ", round(x@dims$resels), "\n\n")
  }
  cat("---- \n\n")
  

  # contrast results----
  cat("## Contrast Results \n\n")
  for (i in seq_len(ncon)) {
    cat("### ", x@C[[i]]$name, "\n\n")
    
    if (class(x@C[[i]]$results) != "character")
      cat(paste("![conimg", i, "](", docname, "_", x@C[[i]]$name, ".png) \n\n", sep = ""))
    
    ## render results----
    cat("Contrast weights: \n\n")
    c <- t(x@C[[i]]$c)
    colnames(c) <- colnames(x@X$X)
    print(c)
    cat("\n\n")
    if (x@control$rft) {
      if (class(x@C[[i]]$results) == "character") {
        cat(x@C[[i]]$results, "\n\n")
      } else {
        cat("#### Set-level \n\n")
        cat("  Clusters = ", ncol(x@C[[i]]$results$clusterLevel), "\n\n")
        cat("  p-value  = ", x@C[[i]]$results$setLevel, "\n\n")
        
        cat("#### Cluster-Level: \n\n")
        print(round(x@C[[i]]$results$clusterLevel, 3))
        cat("\n\n")
        
        cat("#### Peak-level: \n\n")
        print(round(x@C[[i]]$results$peakLevel, 3))
        cat("\n\n")
      }
    } else {
      if (class(x@C[[i]]$results) == "character") {
        cat(x@C[[i]]$results, "\n\n")
      } else {
        print(x@C[[i]]$results)
      }
    }
    cat("Interest degrees of freedom = ", x@C[[i]]$dims$idf, "\n\n")
    cat("      Statistical threshold = ", round(x@C[[i]]$sthresh, 2), "\n\n")
    cat("          Cluster threshold = ", x@C[[i]]$cthresh, "\n\n")
    cat("             Threshold type = ", x@C[[i]]$threshType, "\n\n")
    cat("---- \n\n")
  }
  sink(type = "message")
  sink()
  
  rmarkdown::render(md)
  rmarkdown::render(md, pdf_document())
  return(TRUE)
  # markdown::markdownToHTML(md, paste(docname, ".html", sep = ""))
  # system(paste("pandoc -s ", paste(docname, ".html", sep = ""), "-o ", paste(docname, ".pdf", sep = ""), sep = ""))
}

#' @export
#' @docType methods
#' @details \strong{getCluster} Retrieve a cluster of a specific contrast.
#' @rdname iModel-methods
getCluster <- function(x, contrast, value) {
  out <- antsImageClone(x@C[[contrast]]$clusterImage)
  if (!missing(value))
    out[out == value]
  return(out)
}

#' @export
#' @docType methods
#' @details \strong{plot} Create plot of variables against specific cluster
#'  within contrast.
#' @rdname iModel-methods
setMethod("plot", "iModel", function(x, contrast, cluster) {
  if (length(contrast) > 1)
    stop("Contrast must be of length 1.")
  if (is.null(x@C[[contrast]]$clusterImage)) {
    warning("No results to plot.")
    return(0)
  } else {
    # create mask for specific cluster
    clustimg <- x@C[[contrast]]$clusterImage
    clustimg[clustimg != cluster] <- 0
    clustimg[clustimg != 0] <- 1
    plot(x@iData, x@y, mask = clustimg)
  }
})

#' Control parameters for RFT based analyses
#' 
#' Auxillary function for controlling \code{\link{iModel}} fitting.
#' 
#' @param cf Critical F-threshold for selecting voxels over which the non-sphericity is estimated (default = \code{0.001}).
#' @param scr Logical. scale residuals? (default = \code{TRUE}).
#' @param sar Number of residual images to sample for estimating the FWHM (default = \code{64}).
#' @param n images in conjunction (default = \code{1}).
#' @param iso logical. should images be assumed to be isotropic? (default = \code{TRUE}).
#' @param os offset weighting for iteratively reweighted least squares (default = \code{3}).
#' @param rft logical. should voxels be estimated in resel space for random field theory analysis (default = \code{TRUE}).
#' @param optimization method "none", "IWLS", or "REML" (default = \code{"none"}).
#' @return 
#' 
#' @export iControl
iControl <- function(cf = 0.05, scr = TRUE, sar = 64, n = 1, iso = TRUE, os = 3, rft = TRUE) {
  list(cf = cf,
       scr = scr,
       sar = sar,
       n = n,
       iso = iso,
       os = os,
       rft = rft,
       opt = "none")
}