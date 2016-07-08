#' Summarizing the F-statistic of rftModel fits
#' 
#' 
#' @param object object of class \code{rftModel}
#' @param cthresh cluter threshold for contrasts
#' @param ... optional arguments passed to update.rftModel (see \code{rftModel-class} for more details)
#' 
#' @return Produces an rftContrast object
#' 
#' @details
#' 
#' @author Zachary P. Christensen
#' 
#' @method anova rftModel
#' @describeIn rftContrast-class
anova.rftModel <- function(rftModel, contrastMatrix, cthresh = 150, isotropic = NULL, ...) {
  ncon <- nrow(contrastMatrix)
  out <- list()
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  for (i in seq_len(ncon)) {
    out[[i]] <- rftContrast(rftModel = rftModel, name = connames[i])
    out[[i]] = out[[i]]$setC(contrastMatrix[i, ])
    out[[i]] = out[[i]]$setX1()
    out[[i]] = out[[i]]$fstat()
    out[[i]] = out[[i]]$solve(cthresh, isotropic, ...)
  }
  return(out)
}



#' 
#' 
#' @field rftModel
#' @field C
#' @field X0
#' @field X1
#' @field trMV
#' @field trMVMV
#' @field idf
#' @field fieldType
#' @field contrastImage
#' @field clusterImage
#' @field setLevel
#' @field clusterLevel
#' @field peakLevel
#' @field st statistical threshold
#' @field ct cluster threshold
#' 
#' @export rftContrast
rftContrast <- 
  setRefClass("rftContrast",
              fields = list(
                rftModel = "rftModel",
                C = "numeric",
                X0 = "ANY",
                X1 = "ANY",
                trMV = "numeric",
                trMVMV = "numeric",
                idf = "numeric",
                fieldType = "character",
                contrastImage = "antsImage",
                clusterImage  = "antsImage",
                setLevel = "numeric",
                clusterLevel = "data.frame",
                peakLevel = "data.frame",
                st = "numeric",
                ct = "numeric"),
              methods = list(
                setC = function(con) {
                  'set contrast'
                  if (length(con) != rftModel$dims$npred)
                    stop("the contrast length must be equal to the number of predictors")
                  tmp <- matrix(con, ncol = 1)
                  rownames(tmp) <- colnames(model.matrix(rftModel))
                  C <<- tmp
                },
                setX1 = function() {
                  x <- rftModel$KWX
                  # instea off util_cukx
                  X1$uKX1 <- tcrossprod(diag(x$d[seq_len(x$rk) ]), x$v[, seq_len(x$rk) ]) %*% C
                  X1$uKX1[X1$uKX1 < x$tol] <- 0
                  
                  # instead of util_cukxp
                  X0$uKX0 <- tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]), x$v[, seq_len(x$rk)])
                  X0$uKX0[X0$uKX0 < x$tol] <- 0
                },
                setTrMV = function() {
                  rk <- rftModel$KWX$rk
                  out <- list()
                  
                  out$X <- rftModel$KWX$u[, seq_len(rk)] %*% X1$uKX1
                  svdx <- svd(t(X1$uKX1))
                  out$v <- svdx$v
                  out$u <- svdx$u
                  out$d <- svdx$d
                  out$tol <- max(dim(out$X)) * max(absx$d) * .Machine$double.eps
                  out$rk <- sum(out$d > out$tol)
                  
                  u <- out$u[, seq_len(rk)]
                  Vu <- rftModel$XV %*% u
                  trMV <<- sum(u %*% Vu)
                  trMVMV <<- norm(crossprod(u, Vu), "f")^2
                  idf <<- (trMV^2) / trMVMV
                },
                fStat = function() {
                  'calculate the f-statistic of the contrast'
                  #B <- as.vector(contrast) * as.matrix(rftMod$beta)
                  #mvm <- (t(B %*% rftMod$X)- colSums(as.matrix(rftMod$imageMatrix)))
                  ## mvm <- colSums(t(colSums(t(crossprod(B, crossprod(X, M)) %*% X %*% B))^2)) / trMV
                  #contrastImage <<- mvm / rftMod$mrss
                  
                  # set hsqr
                  out <- list()
                  out$X <- X1$uKX1
                  svdx <- svd(t(X1$uKX1))
                  out$v <<- svdx$v
                  out$u <<- svdx$u
                  out$d <<- svdx$d
                  out$tol <<- max(dim(out$X)) * max(absx$d) * .Machine$double.eps
                  out$rk <<- sum(out$d > out$tol)
                  x <- rftModel$KWX
                  cukx <- tcrossprod(diag(x$d[seq_len(x$rk) ]), x$v[, seq_len(x$rk) ])
                  out <- crossprod(out$u[, seq_len(out$rk)], cukx)
                  
                  # create ESS as sum of squared weighted beta images
                  ess <- colSums((out %*% rftModel$beta)^2) / trMV
                  
                  y <- rftModel$y
                  mask <- rftModel$imgData$imgList[y]$mask
                  contrastImage <<- makeImage(mask, ess / rftModel$mrss)
                  fieldType <<- "F"
                },
                tStat = function() {
                  'calculate the t-statistic of the contrast'
                  Vc <- crossprod(C, rftModel$covBeta) %*% C
                  se <- sqrt(rftMod$mrss * as.numeric(Vc))
                  tvec <- crossprod(C, rftMod$beta) / se
                  contrastImage <<- makeImage(rftMod$mask, tvec)
                  fieldType <<- "T"
                },
                solve = function(cthresh, isotropic = NULL, ...) {
                  'threshold contrastImage and solve for set-level, cluster-level, and peak-level statistics'
                  if (!is.null(isotropic))
                    isotropic <- rftModel$rpvImage
                  else
                    rpvimg <- NULL
                  
                  ll <- c(threshType = "pRFT", pval = .05, pp = .001, n = 1, statdir = NULL, verbose = FALSE)
                  ll[names(...)] <- list(...)
                  
                  z <- rftResults(contrastImage, rftModels$dims$resels, rftModel$dims$fwhm,
                                  c(idf, rftModel$dims$edf), fieldType, rpvImage = isotropic,
                                  k = cthresh, ll$threshType, ll$pval, ll$pp, ll$n, ll$statdir, ll$verbose)
                  
                  threshold <<- z$threshold
                  clusterImage <<- z$clusterImage
                  setLevel <<- z$setLevel
                  peakLevel <<- z$peakLevel
                  clusterLevel <<- z$clusterLevel
                },
                show = function() {
                  'concise print out of rftContrast'
                  cat("rftContrast object: ", name)
                  
                  cat("\nContrast weights: \n")
                  print(t(C))
                  
                  cat("Results \n")
                  cat(" Set-level: \n")
                  cat(" Clusters = ", ncol(clusterLevel), "\n")
                  cat(" p-value = ", setLevel, "\n\n")
                  
                  cat("  Cluster-Level: \n")
                  print(round(clusterLevel, 3))
                  cat("\n")
                  
                  cat(" Peak-level: \n")
                  print(round(peakLevel, 3))
                  cat("\n___\n")
                  
                  cat("Interest degrees of freedom = ", idf, "\n")
                  cat("Residual degrees of freedom = ", rftModel$dims$rdf, "\n")
                  cat("Statistical threshold = ", round(st, 2), "\n")
                  cat("Cluster threshold = ", ct, "\n")
                  cat("Statistical threshold = ", st, "\n")
                  cat("FWHM = ", round(rftModel$dims$fwhm, 2), "\n")
                  cat("Resels = ", round(rftModel$dims$resels), "\n")
                  cat("Voxels = ", rftModel$dims$nvox, "\n")
                }
              )
            )

#' @param object object of class \code{rftModel}
#' @param cthresh cluter threshold for contrasts
#' @param ... optional arguments passed to rftModelUpdate (see \code{rftModelUpdate} for more details)
#' 
#' @return
#' \item{call} {original call to \code{rftFit}}
#' \item{contrastResults} {Contains the resulting image for each contrast along with all the results produce by \code{rftResults} for that particular contrast.}
#' \item{fwhm} {full width at half maxima}
#' \item{resels} {resolution elements}
#' \item{voxels} {number of active voxels (number of voxels in mask)}
#' 
#' @details Summarizes the t-statistic of rftModel fits
#' 
#' @author Zachary P. Christensen
#' 
#' @method summary rftModel
#' @describeIn rftContrast-class
summary.rftModel <- function(rftModel, contrastMatrix, cthresh = 150, isotropic = NULL, ...) {
  ncon <- nrow(contrastMatrix)
  out <- list()
  connames <- rownames(contrastMatrix)
  if (is.null(connames))
    connames <- paste("Contrast_", seq_len(ncon), sep = "")
  for (i in seq_len(ncon)) {
    out[[i]] <- rftContrast(rftModel = rftModel, name = connames[i])
    out[[i]] = out[[i]]$setC(contrastMatrix[i, ])
    out[[i]] = out[[i]]$setX1()
    out[[i]] = out[[i]]$tstat()
    out[[i]] = out[[i]]$solve(cthresh, isotropic, ...)
  }
  return(out)
}
