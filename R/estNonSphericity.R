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
estNonSphericity <- function(x, Y, X, W = NULL, S, mask, df, alpha) {
  
  # add variance component for each level of factor i
  xVi <- x$xVi

  # compute Hsqr and F-threshold under i.i.d.
  if (!is.null(xVi$Fcontrast))
    con <- util_setCon(x, "usc", "F", "c", xVi$Fcontrast, x$KWX)
  else {
    iX0 <- c(x$iB, x$iG)
    con <- util_setCon(x, "eoi", "F", "iX0", iX0, x$KWX)
  }

  if (!is.null(con$c)) {
    X1 <- util_X1(con, x$KWX)
    hsqr <- util_hsqr(con, x$KWX)
    trMV <- util_trMV(X1, oneout = TRUE)
  } else {
    trMv <- 1
    hsqr <- Inf
  }

  # threshold for voxels entering non-sphericity estimates
  ctrl <- x$control
  xVi$UFp <- ctrl$cF
  UF <- qf(1 - ctrl$cF, x$dims$idf, x$dims$rdf)

  good <- (colSums((hsqr %*% beta)^2) / trMV ) > (UF * mrss)
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
      Xp <- lappend(Xp, x$K[[i]]$X0)

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
