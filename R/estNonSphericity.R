## @param imat image matrix
## @param X structure containing design matrix information
## @param W optional whitening/weighting matrix
## @param S structure describing intrinsic non-sphericity
## @param mask structure containing masking information
## @param alpha level for initial threshold
##
## @return
## V estimated non-sphericity, trace(V) = rankd(V)
## h hyperparameters xVi.V = xVi.h(1)*xVi.Vi{1} + ...
## xVi$Cy spatially whitened <Y*Y'> (used by ReML to estimate h)
##
estNonSphericity <- function(object, x, Y, X, W = NULL, S, mask, df, alpha) {
  
  X <- object$X
  if (is.null(X$W))
    X$W <- diag(nrow(X$X))
  xVi <- object$xVi

  # compute Hsqr and F-threshold under i.i.d.
  X$KWX <- .setx(rftFilter(object$K, X$W %*% X$X))
  X$XX <- .pinvx(X$KWX)
  
  if (!is.null(xVi$Fcontrast))
    con <- .setcon("usc", "F", "c", xVi$Fcontrast, X$KWX)
  else 
    con <- .setcon("eoi", "F", "iX0", c(X$iB, X$iG), X$KWX)
  
  if ((!is.null(con[[1]]$c)) | length(con[[1]]$c) == 0) {
    X1 <- .X1(con[[1]], X$KWX)
    hsqr <- .hsqr(con[[1]], X$KWX)
    trmv <- .trMV(X1, oneout = TRUE)
  } else {
    trmv <- 1
    hsqr <- Inf
  }
  trrv <- .trRV(X$KWX, oneout = TRUE)
  
  # threshold for voxels entering non-sphericity estimates
  ctrl <- object$control
  xVi$UFp <- ctrl$cF
  UF <- qf(1 - ctrl$cF, trmv, trrv)

  good <- (colSums((hsqr %*% beta)^2) / trmv) > (UF * mrss)
  q <- sqrt(1 / mrss)[good]
  q <- t(t(imgData[y]$imageMatrix[, good]) * q)
  Cy <- tcrossprod(q)

  # ReML Estimation
  if (is.list(K)) {
    m <- length(xVi$Vi)
    h <- rep(0, m)
    V <- matrix(0, n, n)
    for (i in 1:length(K)) {
      # extract blocks from bases
      q <- K[[i]]$row
      p <- c()
      QP <- list()
      for (j in 1:m) {
        if (any(xVi$Vi[[j]][q, q] != 0)) {
          Qp <- lappend(Qp, xVi$Vi[[j]][q, q])
          p <- c(p, j)
        }
      }

      # design space for ReML (with confounds in filter)
      Xp <- X[q, ]
      Xp <- lappend(Xp, object$K[[i]]$X0)

      # ReML
      reml <- rftModelOptimize(Cy[q, q], Xp, Qp)
      V[q, q] <- V[q, q] + reml$Vp
      h[p] <- reml$hp
    }
  } else {
    reml <- rftModelOptimize(Cy, x$X, xVi$Vi)
    V <- reml$V
    h <- reml$h
  }

  V <- V * n / sum(diag(V))
  xVi$h <- h
  xVi$V <- V
  xVi$Cy <- Cy
  return(xVi)
}
