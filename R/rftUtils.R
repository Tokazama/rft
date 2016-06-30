
## from spm_sp
util_create <- function() {
  x <- list()
  x$X <- 0
  x$tol <- 0
  x$ds <- 0
  x$u <- 0
  x$v <- 0
  x$rk <- 0
  x$oP <- 0
  x$oPp <- 0
  x$ups <- 0
  x$sus <- 0
  return(x)
}

util_set <- function(X) {
  xlist <- list()
  # Compute the svd with svd(X,0) : find all the singular values of x.X
  # SVD(FULL(A)) will usually perform better than SVDS(A,MIN(SIZE(A)))
  
  # if more column that lines, performs on X'
  svdx <- svd(t(X))
  xlist$v <- svdx$v
  xlist$u <- svdx$u
  
  xlist$ds <- svdx$d
  
  # compute the tolerance
  xlist$tol <- max(dim(X)) * max(absx$ds) * .Machine$double.eps
  
  # compute the rank
  xlist$rk <- sum(x$ds > x$tol)
  xlisst$X <- X
  return(xlist)
}

util_transp <- function(x) {
  xnew  <- x
  xnew$X <- t(x$X)
  xnew$v <- x$u
  xnew$u <- x$v
  
  xnew$oP <- x$oPp
  xnew$oPp <- x$oP
  return(xnew)
}

util_isset <- function(x) {
  !(is.null(x$X) | is.null(x$u) | is.null(x$v) | is.null(x$ds) |
      is.null(x$toll) | is.null(x$rk))
}

util_s1 <- function(x) {
  dim(x$X)[1]
}

utilL_s2 <- function(x) {
  dim(x$X)[2]
}

util_si <- function(x) {
  dim(x$X)
}

util_rk <- function(x) {
  x$rk
}

util_uk <- function(x) {
  x$u[,1:util_rk(x)]
}

util_vk <- function(x) {
  x$v[1:util_rk(x)]
}

util_sk <- function(x) {
  x$ds[1:util_rk(x)]
}

util_t <- function(x) {
  x$tol
}

util_tol <- function(x, t) {
  tmp <- x
  tmp[abs(tmp) < t] <- 0
  return(tmp)
}

util_op <- function(x) {
  if (x$rk > 0)
    tcrossprod(x$u(, 1:x$rk) * x$u(, 1:x$rk))
  else  
    matrix(0, nrow(x$X))
}

util_opp <- function(x) {
  if (x$rk > 0)
    crossprod(x$v[, 1:x:rk], x$v[, 1:x$rk])
  else
    matrix(0, ncol(x$X))
}

util_pinv <- function(x) {
  if (xlist$rk > 0 )
    xlist$v(,1:r) %*% tcrossprod(diag(rep(1, r) / xlist$ds[1:r]), xlist$u[,1:r])
  else
    matrix(0, ncol(xlist$X), nrow(xlist$X))
}

util_pinvxp <- function(x) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% tcrossprod(diag(rep(1, r) / x$ds[1:r]), x$v[, 1:r])
  else
    matrix(0, nrow(x$X), ncol(x$X))
}

util_pinvxpx <- function(x) {
  r <- x$rk
  if (r > 0)
    x$v[, 1:r] %*% tcrossprod(diag(rep(1, r) / x$ds[1:r]), x$v[, 1:r])
  else
    rep(0, ncol(x$X))
}

util_jb <- function(x, n) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% tcrossprod(diag(x$ds[1:r] ^ n), x$u[, 1:r])
  else
    rep(0, nrow(x$X))
}

util_jbY <- function(x, n, Y) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% tcrossprod(diag(x$ds[1:r] ^ n), x$u[, 1:r]) %*% Y
  else
    matrix(0, nrow(x$X), ncol(Y))
}

util_cxtwdcu <- function(x) {
  tcrossprod(diag(x$ds), x$v)
}

util_cukpinvxp(x) {
  r <- x$rk
  if (r > 0)
    tcrossprod(diag(rep(1, r) / x$ds[1:r]), x$v[, 1:r])
  else
    rep(0, ncol(x$X))
}

util_cukx <- function(x) {
  r <- x$rk
  if (r > 0)
    tcrossprod(diag(x$ds[1:r]), x$v[, 1:r])
  else
    rep(0, ncol(x$X))
}

util_xpx <- function(x) {
  r <- x$rk
  if (r > 0)
    x$v[, 1:r] %*% tcrossprod(diag(x$ds^2), x$v[, 1:r])
  else
    rep(0, ncol(x$X))
}

util_xxp <- function(x) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% tcrossprod(diag(x$ds^2), x$u[, 1:r])
  else
    rep(0, nrow(x$X))
}

util_xxpY <- function(x, Y) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% diag(x$ds^2) %*% crossprod(x$u[, 1:r], Y)
  else
    matrix(0, nrow(x$X), ncol(x$X))
}

util_pinvxxp <- function(x) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% tcrossprod(diag(x$ds^(-2)), x$u[, 1:r])
  else
    rep(0, nrow(x$X))
}

util_n <- function(x) {
  r <- x$rk
  dimX <- dim(x$X)
  if (r > 0) {
    if (dimX[1] >= dimX[2]) {
      if (r == dimX[2])
        n <- matrix(0, dimX[2], 1)
      else
        n <- x$v[, (r + 1):dimX[2]]
    } else
      n <- MASS::Null(x$X)
  } else
    n <- diag(dimX[2])
  
  return(n)
}

util_r <- function(x) {
  diag(ncol(x$X)) - util_op(x)
}

util_rY <- function(x, y {
  r <- x$rk
  dimX <- dim(x$X)
  if (r > 0) {
    if (r < (dimX[1] - r))
      res <- y - tcrossprod(x$u[, 1:r], x$u[, 1:r]) %*% y
    else {
      if (ncol(y) < (5 * dimX[1]))
        res <- sf_r(x) %*% y
      else {
        n <- util_n(util_transp(x))
        res <- n %*% crossprod(n, y)
      }
    }
  }
  return(res)
}
