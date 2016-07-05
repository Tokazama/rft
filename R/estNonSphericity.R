# @param imat image matrix
# @param X structure containing design matrix information
# @param W optional whitening/weighting matrix
# @param S structure describing intrinsic non-sphericity
# @param mask structure containing masking information
# @param alpha level for initial threshold
# 
# @return 
# V estimated non-sphericity, trace(V) = rankd(V)
# h hyperparameters xVi.V = xVi.h(1)*xVi.Vi{1} + ...
# xVi.Cy spatially whitened <Y*Y'> (used by ReML to estimate h)
#
estNonSphericity <- function(x, Y, X, W = NULL, S, mask, df, alpha) {
  
  # add variance component for each level of factor i
  
  n <- dim(x$dm)[1]
  rank <- dim(X)[2]
  nvox <- dim(Y)[2]
  # compute Hsqr and F-threshold under i.i.d.
  
  # threshold for voxels--------------------------------------------------------
  UF <- qf(1 - alpha, x$df[1], x$df[2])
  
  # remove filter confounds-----------------------------------------------------
  # KWY <- frequencyFilterfMRI(boldmat, tr, freqLo = 0.01, freqHi = 0.1, opt = "butt")
  KWY <- Y * x$W
  
  # OLS estimation--------------------------------------------------------------
  z <- .lm.fit(x$dm, KWY)
  rss <- rowSums(z$residuals^2)
  
  # F-threshold and acccumulate spatially whitened Y*Y'
  se <- t(as.matrix(sqrt((rss / x$dof[2]) * (contrast %*% object$XX %*% contrast))))
  se[se == 0] <- 1 # control for NULLS that result from dividing by 0
  j <- ((contrast %*% z$coefficients) / se)^2 > (UF * rss / trRV) # logical vector of voxels above threshold
  rm(z, se)
  
  
  sum((Hsqr * beta).^2, 1)/trMV > UF * ReSS / trRV
  
  # Pooled Variance Estimation--------------------------------------------------
  voxseq <- seq(1:rftmod$volumes$voxels, by = controlvals$chunkSize)
  voxseq_iter <- length(voxseq)
  s <- 0 # number of groups with voxels above threshold
  for (i in 1:nchunks) {
    if (i == nchunks)
      subvox <- voxseq[i]:rftmod$volumes$voxels
    else
      subvox <- voxseq[i]:(voxseq[i + 1] - 1)
    
    tmpy <- object$y[, subvox]
    
    # voxels above f-threshold
    olsfit <- .lm.fit(rftmod$parameters$X, tmpy)
    rss <- colSums(olsfit$residuals^2)
    fstat <- (((olsfit$fitted.values - colMeans(tmpy))^2) / object$dof[1]) / (rss / object$dof[2])
    good <- fstat > fthresh
    
    if (any(good != FALSE)) {
      q <- diag(sqrt(trRV / rss[, good]))
      tmpy  <- tmpy %*% q
      Cy <- Cy + tcrossprod(tmpy)
      s <- s + 1
    }
  }
  if (s < 1)
    stop("There are no significant voxels")
  Cy <- Cy / s
  
  # ReML Estimation-------------------------------------------------------------
  if (is.list(K)) {
    m <- length(Vi)
    h <- rep(0, m)
    V <- matrix(0, n, n)
    for (i in 1:length(K)) {
      # extract blocks from bases
      q <- nrow(K[[i]]$row)
      p <- 0
      QP <- list()
      for (j in 1:m) {
        if (nnz(xVi.Vi[[j]][q, q]) {
          Qp <- lappend(Qp, Vi[[j]][q, q]);
          p <- c(p, j);
        }
      }
      
      # design space for ReML (with confounds in filter)
      Xp <- X[q, ]
      Xp = c(Xp, xX.K[[i]].X0)
      
      # ReML
      reml <- rftRML(Cy(q, q), Xp, Qp)[c("Vp", "hp")]
      V(q, q) <- V(q, q) + reml$Vp
      h(p) <- reml$hp
      
    }
  } else
    reml = rftRML(Cy, X, Vi)[c("V", "h")]
  
  V = reml$V * n /sum(diag(reml$V))
  xVi.h = reml$h
  xVi.v = v
  xVi.Cy = Cy
}
