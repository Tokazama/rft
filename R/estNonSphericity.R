# Estimate non-sphericity of iModel
# 
# Mainly for internal use
# 
# @param object 
# 
# @return changes the following values in the \code{xVi} slot of the iModel object
# \item{V} {estimated non-sphericity, trace(V) = rankd(V)}
# \item{h} {hyperparameters xVi.V = xVi.h(1)*xVi.Vi{1} + ...}
# \item{Cy} {spatially whitened <Y*Y'> (used by ReML to estimate h)}
# 
# @export estNonSphericity
estNonSphericity <- function(object) {
  if (!is.null(object@xVi$Fcontrast))
    con <- .setcon("usc", "F", "c", object@xVi$Fcontrast, object@X$KWX)
  else 
    con <- .setcon("eoi", "F", "iX0", c(object@X$iB, X$iG), object@X$KWX)
  
  if ((!is.null(con[[1]]$c)) | length(con[[1]]$c) == 0) {
    X1 <- .X1(con[[1]], object@X$KWX)
    hsqr <- .hsqr(con[[1]], object@X$KWX)
    trmv <- .trMV(X1, oneout = TRUE)
  } else {
    trmv <- 1
    hsqr <- Inf
  }
  trrv <- .trRV(object@X$KWX, oneout = TRUE)
  
  # threshold for voxels entering non-sphericity estimates
  object@xVi$UFp <- object@control$cf
  UF <- qf(1 - object@control$cf, trmv, trrv)
  
  
  end <- object@control$chunksize
  nchunks <- floor(fit@dims$nvox / end)
  start <- 1 + end
  Cy <- matrix(0, object@dims$nimg, object@dims$nimg)
  s <- 0
  for (j in seq_len(nchunks)) {
    if (j == nchunks)
      vrange <- start:fit@dims$nvox
    else
      vrange <- start:end
    
    KWY <- iFilter(object@X$K, object@X$W %*% object@iData[[object@y]]@iMatrix[, vrange])
    object@beta[, vrange] <- object@X$XX %*% object@KWY
    object@res[, vrange] <- .res(X$KWX, KWY)
    object@mrss[, vrange] <- colSums(object@res[, vrange]^2) / object@X$trRV
    if (sr)
      object@res[, vrange] <- t(t(object@res[, vrange]) * (1 / as.numeric(object@mrss[, vrange])))
    
    start <- start + object@control$chunksize
    end <- end + object@control$chunksize
    
    good <- vrange[(colSums((hsqr %*% object@beta[, vrange])^2) / trmv) > (UF * object@mrss[, vrange])]

    if (length(good) > 0) {
      q <- sqrt(1 / mrss[, good])
      q <- t(t(object@iData[[object@y]]@imageMatrix[, good]) * q)
      Cy <- tcrossprod(q)
      s <- s + length(good)
    }
  }
  if (s > 0)
    Cy <- Cy / s
  else
    warning("No voxels appear to be significant.")
  
  if (is.list(object@X$K)) {
    m <- length(object@xVi$Vi)
    h <- rep(0, m)
    V <- matrix(0, n, n)
    for (i in seq_len(length(object@X$K))) {
      # extract blocks from bases
      q <- object@X$K[[i]]$row
      p <- c()
      QP <- list()
      for (j in seq_len(m)) {
        if (any(object@xVi$Vi[[j]][q, q] != 0)) {
          Qp <- lappend(Qp, object@xVi$Vi[[j]][q, q])
          p <- c(p, j)
        }
      }
      
      # design space for ReML (with confounds in filter)
      Xp <- X[q, ]
      Xp <- lappend(Xp, object@X$K[[i]]$X0)
      
      # ReML
      reml <- iREML(Cy[q, q], Xp, Qp, mi = object$control$mi)
      V[q, q] <- V[q, q] + reml$Vp
      h[p] <- reml$hp
    }
  } else {
    reml <- iREML(Cy, object@X$X, object@xVi$Vi, mi = object$control$mi)
    V <- reml$V
    h <- reml$h
  }
  
  V <- V * n / sum(diag(V))
  
  object@xVi$h <- h
  object@xVi$V <- V
  object@xVi$Cy <- Cy
  return(x)
}