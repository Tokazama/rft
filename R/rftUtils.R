## from spm_sp
## to do:
## isspc
## issetspc
## size = just gives the dimensions of x$X
isinsp = function(x, c, tol) {
  'is in space or is in dual space'
  if (missing(tol))
    tol <- x$tol
  if (nrow(x$X) != nrow(c)) {
    warning('Vector dim don''t match col dim : not in space')
    return(0)
  } else {
    out <- all(abs(util_op(x) %*% c - c) <= tol)
    return(colSums(out) > 0)
  }
}

isinspp = function(x, c, tol) {
  'is in space or is in dual space'
  if (missing(tol))
    tol <- x$tol
  if (ncol(x$X) != nrow(c)) {
    warning('Vector dim don''t match row dim : not in space')
    return(0)
  } else {
    out <- all(abs(util_opp(x) %*% c - c) <= tol)
    return(colSums(out) > 0)
  }
}

eachinsp = function(x, c, tol) {
  'each column of c in space or in dual space'
  if (missing(tol))
    tol <- x$tol
  if (nrow(x$X) != nrow(c)) {
    warning('Vector dim don''t match col. dim : not in space')
    return(0)
  } else
    return(all(abs(util_op(x) %*% c - c) <= tol))
}

eachinspp = function(x, c, tol) {
  'each column of c in space or in dual space'
  if (missing(tol))
    tol <- x$tol
  if (ncol(x$X) != nrow(c)) {
    warning('Vector dim don''t match row. dim : not in space')
    return(0)
  } else
    return(all(abs(util_opp(x) %*% c - c) <= tol))
}

util_equal = function(x, X2) {
  'test wether two spaces are the same'
  x2 <- util_set(X2)
  maxtol <- max(x$tol, x2$tol)
  return(all(util_isinsp, x, X2, maxtol) && all(util_isinsp(x2, x$X, maxtol))
}

check = function(x, y, tol) {
  'set y and tol according to arguments'
  if (missing(tol))
    tol <- x$tol
  if (missing(y))
    return(util_tol(x, tol))
  else
    return(util_tol(y, tol))
}

util_nop = function(x, y, check) {
  'project into null space'
  n <- util_n(x)
  if (missing(y))
    n <- tcrossprod(n)
  else
    n <- n %*% crossprod(n, y)
  if (check)
    return(util_tol(n, x$tol))
  else
    return(n)
}

util_res = function(x, y, check) {
  if (missing(y))
    out <- util_r(x)
  else
    out <- util_rY(x, y)
  if (check)
    return(util_tol(out, x$tol))
  else
    return(out)
}

util_carrotp <- function(x, n, y, check) {
  if (missing(n))
    n <- 1
  if (missing(y))
    n <- util_jbp(x, n)
  else
    n <- util_jbp(x, n) %*% y
  if (check)
    return(util_tol(n, x$tol))
  else
    return(n)
}

util_carrot <- function(x, n, y, check) {
  if (missing(n))
    n <- 1
  if (missing(y))
    out <- util_jb(x, n)
  else
    out <- util_jb(x, n, y)
  if (check)
    return(util_tol(out, x$tol))
  else
    return(out)
}

util_np <- function(x) {
  'null space of x transposed'
  util_n(util_transp(x))j
}

util_opColon <- function(x, y) {
  if (missing(y))
    util_tol(util_op(x), x$tol)
  else
    util_tol(util_op(x) %*% y, x$tol)
}

util_oppColon <- function(x, y) {
  if (missing(y))
    util_tol(util_opp(x), x$tol)
  else
    util_tol(util_opp(x) %*% y, x$tol)
}

util_xDash <- function(x, y) {
  if (missing(y))
    util_pinv(x)
  else
    util_pinv(x) %*% y
}

util_xDashColon <- function(x, y) {
  if (missing(y))
    util_tol(util_pinv(x), util_t(x))
  else
    util_tol(util_pinv(x) %*% y, util_t(x))
}

util_xpDash <- function(x, y) {
  if (missing(y))
    util_pinvxp(x)
  else
    util_pinvxp(x) %*% y
}

util_xpDashColon <- function(x, y) {
  if (missing(y))
    util_tol(util_pinvxp(x), util_t(x))
  else
    util_tol(util_pinvxp(x) %*% y, util_t(x))
}

util_cukxpDash <- function(x, y) {
  if (missing(y))
    util_cukpinvxp(x)
  else
    util_cukpinvxp(x) %*% y
}

util_cukxpDashColon <- function(x, y) {
  if (missing(y))
    util_tol(util_cukpinvxp(x), x$tol)
  else
    util_tol(util_cukpinvxp(x) %*% y, x$tol)
}

util_cukxColon <- function(x, y) {
  if (missing(y))
    util_tol(util_cukx(x), x$tol)
  else
    util_tol(util_cukx(x, y), x$tol)
}

util_ox <- function(x) {
  if (util_rk(x) > 0)
    util_uk(x)
  else
    matrix(0, util_s1(x), 1)
}

util_oxp <- function(x) {
  if (util_rk(x) > )
    util_vk(x)
  else
    matrix(0, util_s2(x), 1)
}

util_xi <- function(x, i) {
  x$X[, i]
}

util_xpi <- function(x, i) {
  t(x$X[i, ])
}

util_cx2cu <- function(x, y) {
  if (missing(y))
    util_cxtwdcu(x)
  else
    util_cxtwdcu(x) %*% y
}

util_cx2cuColon <- function(x, y) {
  if (missing(y))
    util_tol(util_cxtwdcu(x), x$tol)
  else
    util_tol(util_cxtwdcu(x) %*% y, x$tol)
}

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

util_uk <- function(x, y) {
  if (missing(y))
    x$u[,1:util_rk(x)]
  else
    x$u[,1:util_rk(x)] %*% y
}

util_ukColon <- function(x, y) {
  if (missing(y))
    util_tol(util_uk(x), x$tol)
  else
    util_tol(util_uk(x, y), x$tol)
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

util_op <- function(x, y) {
  if (x$rk > 0)
    out <- tcrossprod(x$u(, 1:x$rk) * x$u(, 1:x$rk))
  else  
    out <- matrix(0, nrow(x$X))
  if (missing(y))
    return(y)
  else
    return(out %*% y)
}

util_opp <- function(x, y) {
  if (x$rk > 0)
    out <- crossprod(x$v[, 1:x:rk], x$v[, 1:x$rk])
  else
    out <- matrix(0, ncol(x$X))
  if (missing(y))
    return(y)
  else
    return(out %*% y)
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

util_pinvxpx <- function(x, y) {
  r <- x$rk
  if (r > 0)
    out <- x$v[, 1:r] %*% tcrossprod(diag(rep(1, r) / x$ds[1:r]), x$v[, 1:r])
  else
    out <- rep(0, ncol(x$X))
  if (missing(y))
    return(out)
  else
    return(out %*% y)
}

util_pinvxpxColon <- function(x, y) {
  if (missing(y))
    util_tol(util_pinvxpx(x), x$tol)
  else
    util_tol(util_pinvxpx(x, y), x$tol)
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

util_cukx <- function(x, y) {
  r <- x$rk
  if (r > 0)
    out <- tcrossprod(diag(x$ds[1:r]), x$v[, 1:r])
  else
    out <- rep(0, ncol(x$X))
  if (missing(y))
    return(out)
  else
    return(out %*% y)
}

util_xpx <- function(x, y) {
  r <- x$rk
  if (r > 0)
    out <- x$v[, 1:r] %*% tcrossprod(diag(x$ds^2), x$v[, 1:r])
  else
    out <- rep(0, ncol(x$X))
  if (missing(y))
    return(out)
  else
    return(out %*% y)
}

util_xpxColon <- function(x, y) {
  if (missing(y))
    util_tol(util_xpx(x), x$tol)
  else
    util_tol(util_xpx(x, y), x$tol)
}

util_xxp <- function(x, y) {
  r <- x$rk
  if (r > 0)
    out <- x$u[, 1:r] %*% tcrossprod(diag(x$ds^2), x$u[, 1:r])
  else
    out <-rep(0, nrow(x$X))
  if (missing(y))
    return(out)
  else
    return(out %*% y)
}

util_xxpColon <- functioni(x, y) {
  if (missing(y))
    util_tol(util_xxp(x), x$tol)
  else
    util_tol(util_xxp(x, y), x$tol)
}

util_xxpY <- function(x, Y) {
  r <- x$rk
  if (r > 0)
    x$u[, 1:r] %*% diag(x$ds^2) %*% crossprod(x$u[, 1:r], Y)
  else
    matrix(0, nrow(x$X), ncol(x$X))
}

util_pinvxxp <- function(x, y) {
  r <- x$rk
  if (r > 0)
    out <- x$u[, 1:r] %*% tcrossprod(diag(x$ds^(-2)), x$u[, 1:r])
  else
    out <- rep(0, nrow(x$X))
  if (missing(y))
    return(out)
  else
    return(out %*% y)
}

util_pinvxxpColon <- function(x, y) {
  if (missing(y))
    util_tol(util_pinvxxp(x), x$tol)
  else
    util_tol(util_pinvxxp(x, y), x$tol)
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

util_jbp = functioni(x, n) {
  r <- x$rk
  if (r > 0)
    x$v[, 1:r] %*% tcrossprod(diag(x$ds[1:r]^(n)), x$v[, 1:r])
  else
    rep(0, ncol(x$X))
}

## from spm_SpUtil
util_X0_2_c <- function(X0, x) {
  if (is.null(X0)) {
    sc0 <- utilSet(utilX(x, X0))
    if (sc0$rk)
      c <- utiloPp(x, utilR(sc0))
    else
      c <- utiloPp(x)
    
    
  }
}
