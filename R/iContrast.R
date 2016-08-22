

#' @export
#' @docType methods
#' @description \strong{summary} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-contrasts
setMethod("anova", "iModel", function(object, contrastMatrix, cthresh = 150, verbose = TRUE) {
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(rep(1, ncol(object@X$X)), nrow = 1)
    colnames(contrastMatrix) <- colnames(object@X$X)
    if (any(colnames(contrastMatrix) == "(Intercept)"))
      contrastMatrix[, "Intercept"] <- 0
  }
  
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
    object@C[[i]] <- .setcon(connames[i], fieldType = "F", action, contrast, object@X$KWX)
    
    # X1 <<- .X1(.self, iModel$X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)
    names(object@C[[i]]$dims) <- c("trMV", "trMVMV", "idf")
    
    # solve contrast
    h <- .hsqr(object@C[[i]], object@X$KWX)
    ss <- (rowSums((h %*% object@beta)^2) / object@C[[i]]$dims$trMV)
    object@C[[i]]$contrastImage <- makeImage(mask, ss / object@mrss)
    
    if (!object@control$rft) {
      if (!object@control$iso)
        isotropic <- object$rpvImage
      else
        isotropoic <- NULL
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm,
                      c(object@C[[i]]$dims$idf, object$dims$edf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = cthresh[i], ll$threshType, ll$pval, ll$pp, ll$n, ll$statdir, verbose = verbose)
      
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
  names(object@C)[start:ncon] <- connames
  return(object)
})

#' @export
#' @docType methods
#' @description \strong{summary} Compute f-statistic or t-statistic contrast for iModel object.
#' @rdname iModel-contrasts
setMethod("summary", "iModel", function(object, contrastMatrix, cthresh = 150, verbose = TRUE) {
  if (missing(contrastMatrix)) {
    contrastMatrix <- matrix(rep(1, ncol(object@X$X)), nrow = 1)
    colnames(contrastMatrix) <- colnames(object@X$X)
    if (any(colnames(contrastMatrix) == "(Intercept)"))
      contrastMatrix[, "Intercept"] <- 0
  }
  
  # check arguments----
  if (ncol(contrastMatrix) != object@dims$npred)
    stop("The contrast length must be equal to the number of predictors")
  ncon <- nrow(contrastMatrix)
  if (length(cthresh) != ncon)
    cthresh <- rep(cthresh, ncon)
  
  # set isotropic option
  if (!object@control$iso)
    isotropic <- object@dims$rpvImage
  else
    isotropic <- NULL
  
  # set contrast names
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  

  mask <- object@iData@iList[[object@y]]@mask
  
  # ensure that if previous contrasts were set they won't be erased
  old <- length(object@C)
  conseq <- seq_len(ncon) + old
  
  ll <- object@control
  colnames(contrastMatrix) <- colnames(object@X$X)
  for (i in conseq) {
    c <- matrix(contrastMatrix[i, ], ncol = 1)
    object@C[[i]] <- .setcon(connames[i - old], fieldType = "T", "c+", c, object@X$KWX)
    
    X1 <- .X1(object@C[[i]], object@X$KWX)
    object@C[[i]]$dims <- .trMV(X1, object@X$V)

    # calculate t-contrast----
    Vc <- crossprod(object@C[[i]]$c, object@X$betaCov) %*% object@C[[i]]$c
    se <- sqrt(object@mrss[] * as.numeric(Vc))
    se[se == 0] <- 1  # prevent formation of NaN
    tvec <- crossprod(object@C[[i]]$c, object@beta[]) / se
    object@C[[i]]$contrastImage <- makeImage(mask, tvec)
    object@C[[i]]$cthresh <- cthresh[i - old]
    
    # summary of information----
    if (object@control$rft) {
      z <- rftResults(object@C[[i]]$contrastImage, object@dims$resels, object@dims$fwhm, 
                      c(object@C[[i]]$dims$idf, object@X$rdf), object@C[[i]]$fieldType, rpvImage = isotropic,
                      k = object@C[[i]]$cthresh, ll$tt, ll$pval, ll$pp, ll$n, NULL, verbose = verbose)
      
      object@C[[i]]$sthresh <- z$threshold
      object@C[[i]]$clusterImage <- z$clusterImage
      object@C[[i]]$results$setLevel <- z$setLevel
      object@C[[i]]$results$peakLevel <- z$peakLevel
      object@C[[i]]$results$clusterLevel <- z$clusterLevel
    } else {
      # don't use random field theory for summary----
      object@C[[i]]$sthresh <-
        statFieldThresh(object@C[[i]]$contrastImage, ll$pval, object@dims$nvox,
                        ll$n, fwhm = c(1, 1, 1), resels(mask, c(1, 1, 1)),
                        c(object@C[[i]]$dims$idf, object$dims$rdf),
                        object@C[[i]]$fieldType, ll$tt, ll$pp)
      
      object@C[[i]]$clusterImage <-
        labelClusters(object@C[[i]]$contrastImage, minClusterSize = object@C[[i]]$cthresh,
                      minThresh = object@C[[i]]$sthresh, maxThresh = Inf)
      
      object@C[[i]]$results <- labelStats(object@C[[i]]$contrastImage, object@C[[i]]$clusterImage)
    }
  }
  names(object@C)[conseq] <- connames
  return(object)
})

# #' @export
# #' @docType methods
# #' @description
# #' @rdname iModel-contrasts
# setMethod("plot", "iModel", function(x, conname) {
#   for (i in seq_len(length(conname))) {
#     plot(x@iData@iList[[x@y]]@mask, x@C[[conname[i]]]$clusterImage)
#   }
# })