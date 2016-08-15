# master doc

#'
#'
#' @slot iModel
#' @slot C
#' @slot X0
#' @slot X1
#' @slot trMV
#' @slot trMVMV
#' @slot idf
#' @slot fieldType
#' @slot contrastImage
#' @slot clusterImage
#' @slot setLevel
#' @slot clusterLevel
#' @slot peakLevel
#' @slot sthresh statistical threshold
#' @slot cthresh cluster threshold
#' 
#' @export iContrast
iContrast <- setClass("iContrast",
                      slots = list(
                        name = "character",
                        iModel = "iModel",
                        c = "matrix",
                        X1 = "ANY",
                        X0 = "ANY",
                        iX0 = "ANY",
                        dims = "list",
                        fieldType = "character",
                        contrastImage = "antsImage",
                        clusterImage  = "antsImage",
                        results = "ANY",
                        sthresh = "numeric",
                        cthresh = "numeric")
                      )

setMethod("initialize", "iContrast", function(.Object, name, fieldType, contrast, model) {
  if (missing(model)) {
    .Object@iModel <- iModel()
    .Object@X1 <- matrix(0, 0, 0)
    .Object@X0 <- matrix(0, 0, 0)
    .Object@c <- matrix(0, 0, 0)
    
    .Object@dims <- list()
    .Object@dims$trMV <- 0
    .Object@dims$trMVMV <- 0
    .Object@dims$idf <- 0
    
    .Object@results <- list()
    .Object@results$clusterLevel <- data.frame()
    .Object@results$peakLevel <- data.frame()
    .Object@results$setLevel <- 0
  } else {
    contrast <- if (missing(contrast)) matrix(0, 0, 0)
    if (length(contrast) != iModel$dims$npred)
      stop("the contrast length must be equal to the number of predictors")
    tmp <- matrix(con, ncol = 1)
    rownames(tmp) <- colnames(iModel@X$X)
    c <- tmp
    xCon <- .setcon(name, fieldType, action, contrast, iModel$X$KWX)
    .Object@name <- xCon$name
    .Object@c <- xCon$c
    .Object@X1 <- xCon$X1
    .Object@X0 <- xCon$X0
    .Object@iX0 <- xCon$iX0
    
    # X1 <<- .X1(.self, iModel$X$KWX)
    out <- .trMV(X1, iModel$X$V)
    .Object@dims <- list()
    .Object@dims$trMV <- out[[1]]
    .Object@dims$trMVMV <- out[[2]]
    .Object@dims$idf <- out[[3]]
  }
  
  # solve contrast
  if (.Object@fieldType == "T") {
    Vc <- crossprod(.Object@c, .Object@iModel@covBeta) %*% .Object@c
    se <- sqrt(.Object@iModel@mrss * as.numeric(Vc))
    tvec <- crossprod(.Object@c, .Object@iModel@beta) / se
    .Object@contrastImage <- makeImage(.Object@iModel@mask, tvec)
  } else if (.Object@fieldType == "F") {
    h <- .hsqr(xCon, .Object@iModel@X$KWX)
    ss <- (rowSums((h %*% .Object@iModel@beta)^2) / .Object@dims$trMV)
    mask <- .Object@iModel@iData@iList[[.Object@iModel@y]]@mask
    .Object@contrastImage <- makeImage(mask, ss / .Object@iModel@mrss)
  }
  
  # get results
  ll <- iModel$control
  .Object@ct <- cthresh
  if (!is.null(!ll$rft)) {
    if (!ll$iso)
      isotropic <- .Object@iModel$rpvImage
    else
      rpvimg <- NULL
    z <- rftResults(.Object@contrastImage, .Object@iModels$dims$resels, .Object@iModel$dims$fwhm,
                    c(.Object@idf, .Object@iModel$dims$edf), .Object@fieldType, rpvImage = isotropic,
                    k = ct, ll$threshType, ll$pval, ll$pp, ll$n, ll$statdir, ll$verbose)
    .Object@sthresh <- z$threshold
    .Object@clusterImage <- z$clusterImage
    .Object@results$setLevel <- z$setLevel
    .Object@results$peakLevel <- z$peakLevel
    .Object@results$clusterLevel <- z$clusterLevel
  } else {
    mask <- .Object@iModel$iData$iList[[iModel$y]]$mask
    .Object@sthresh <- statFieldThresh(contrastImage, ll$pval, .Object@iModel$dims$nvox, ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)), c(.Object@dims$idf, .Object@iModel$dims$rdf), .Object@fieldType, ll$tt, ll$pp)
    
    .Object@clusterImage <- labelClusters(.Object@contrastImage, minClusterSize = .Object@cthresh, minThresh = .Object@sthresh, maxThresh = Inf)
    .Object@results <- labelStats(.Object@contrastImage, .Object@clusterImage)
  }
  
  return(.Object)
})

setMethod("show", "iContrast", function(object) {
  'concise print out of iContrast'
  cat("iContrast object: ", object@name)
  
  cat("\nContrast weights: \n")
  print(t(object@c))
  
  cat("Results \n")
  if (object@iModel$control$rft) {
    cat(" Set-level: \n")
    cat(" Clusters = ", ncol(object@results$clusterLevel), "\n")
    cat(" p-value = ", object@results$setLevel, "\n\n")
    
    cat("  Cluster-Level: \n")
    print(round(object@results$clusterLevel, 3))
    cat("\n")
    
    cat(" Peak-level: \n")
    print(round(object@results$peakLevel, 3))
  } else {
    print(results)
  }
  cat("\n___\n")
  
  cat("Interest degrees of freedom = ", object@dims$idf, "\n")
  cat("Residual degrees of freedom = ", object@iModel$dims$rdf, "\n")
  cat("Statistical threshold = ", round(object@st, 2), "\n")
  cat("Cluster threshold = ", object@cthresh, "\n")
  cat("Statistical threshold = ", sthresh, "\n")
  if (object@iModel$control$rft && length(object@iModel$dims$fwhm > 0)) {
    cat("FWHM = ", round(object@iModel$dims$fwhm, 2), "\n")
    cat("Resels = ", round(object@iModel$dims$resels), "\n")
  }
  cat("Voxels = ", object@iModel$dims$nvox, "\n")
})


#' @details Compute the analysis of variance f-statistic for an iModel object
#' @describeIn iContrast
setMethod("anova", "iModel", function(object, contrastMatrix, cthresh = 150) {
  ncon <- nrow(contrastMatrix)
  out <- list()
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  for (i in seq_len(ncon))
    out[[i]] <- iContrast(connames[i], fieldType = "F", matrix(contrastMatrix[i, ], ncol = 1), object)
  return(out)
})

#' @details Compute the t-statistic for an iModel object
#' @describeIn iContrast
setMethod("summary", "iModel", function(object, contrastMatrix, cthresh = 150) {
  ncon <- nrow(contrastMatrix)
  out <- list()
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  for (i in seq_len(ncon))
    out[[i]] <- iContrast(connames[i], fieldType = "T", matrix(contrastMatrix[i, ], ncol = 1), object)
  return(out)
})