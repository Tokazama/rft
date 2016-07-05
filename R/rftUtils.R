## from spm_sp
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

util_xi <- function(x, i) {
  x$X[, i]
}

util_xpi <- function(x, i) {
  t(x$X[i, ])
}

util_isset <- function(x) {
  !(is.null(x$X) | is.null(x$u) | is.null(x$v) | is.null(x$ds) |
      is.null(x$toll) | is.null(x$rk))
}

util_c2tsp = function(x, c, oneout = FALSE, plus = FALSE) {
  'orthonormal partitioning implied by F-contrast'
  if (missing(x))
    x <- KWX
  if (oneout) {
    if (is.list(c) && is.null(c)) {
      if (plus)
        out <- .self$getCUKXt(x, c, check = FALSE)
      else
        out <- x$X %*% c
    } else if (!is.list(c)) {
      if (plust)
        out <- .self$getCUKXp(x, c, check = TRUE)
      else
        out <- .self$getXp(x, c, check = TRUE)
    }
  } else {
    if (is.list(c) && is.null(c)) {
      if (plus)
        out <- c(.self$getCUKXt(x, c), .self$getCUKX(x, .self$getR(util_set(c))))
      else {
        out <- .self$getXp(x, c, check = FALSE)
        out[out < x$tol] <- 0
        out <- c(out, x$X %*% .self$getR(.self$setKWX(c)))
      }
    } else if (!is.list(c)) {
      if (plust) {
        out <- .self$getCUKX(x, c)
        out[out < x$tol] <- 0
        out2 <- .self$getCUKX(x)
        out2[out2 < x$tol] <- 0
        out <- c(out, out2)
      } else
        out <- c(x$X %*% c, x$X)
    }
  }
  return(out)
}

util_cX1 = function(x, c, fieldType) {
  if (missing(x))
    x <- KWX
  
  out <- list()
  out$name <- name
  out$fieldType <- fieldType
  out$c <- matrix(c, ncol = 1)
  
  tmp <- util_c2tsp(x, c, plus = TRUE)
  out$X1$uKX1 <- tmp[[1]]
  out$X0$uKX0 <- tmp[[2]]
  
  X1 <- util_ox() %*% out$X1$uKX1
}

'coordinates in the basis of X to basis u'
util_cx2cu <- function(x, y, check = FALSE) {
  if (mussing(y))
    out <- tcrossprod(diag(x$d), x$v)
  else
    out <- tcrossprod(diag(x$d), x$v) %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'null space of t(X)'
util_np = function(x) {
  out  <- list()
  out$X <- t(x$X)
  out$v <- x$u
  out$u <- x$v
  
  out$oP <- x$opp
  out$oPp <- x$op
  
  dimX <- dim(out$X)
  if (x$rk > 0) {
    if (dimX[1] >= dimX[2]) {
      if (out$rk == dimX[2])
        n <- matrix(0, dimX[2], 1)
      else
        n <- out$v[, (x$rk):dimX[2]]
    } else
      n <- MASS::Null(out$X)
  } else
    n <- diag(dimX[2])
  return(n)
}

'coordinates of pseudo-inverse of t(X) in the base of uk'
util_cukxp = function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'coordinates of X in the base of uk'
util_cukx = function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- tcrossprod(diag(x$d[seq_len(x$rk) ]), x$v[, seq_len(x$rk) ])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[out < x$tol] <- 0
  return(out)
}

'project into null space X'
util_nop = function(x, y, check = FALSE) {
  dimX <- dim(x$X)
  if (x$rk > 0) {
    if (dimX[1] >= dimX[2]) {
      if (x$rk == dimX[2])
        n <- matrix(0, dimX[2], 1)
      else
        n <- x$v[, (x$rk):dimX[2]]
    } else
      n <- MASS::Null(x$X)
  } else
    n <- diag(dimX[2])
  if (missing(y))
    out <- tcrossprod(n)
  else
    out <- n %*% crossprod(n, y)
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

'projection into null space of t(X)'
util_nopp = function(x, y, check = FALSE) {
  out <- list()
  out$X <- t(x$X)
  out$v <- x$u
  out$u <- x$v
  out$tol <- x$tol
  out$rk <- x$rk
  if (missing(y))
    .self$getNOP(out, check = check)
  else
    .self$getNOP(out, y, check = check)
}

'orthogonal projectors on space of X'
util_op = function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- tcrossprod(x$u(, seq_len(x$rk)), x$u(, seq_len(x$rk)))
  else  
    out <- matrix(0, nrow(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

'orthogonal projectors on space of t(X)'
util_opp = function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- crossprod(x$v[, seq_len(x$rk)], x$v[, seq_len(x$rk)])
  else
    out <- matrix(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

'orhtonormal basis sets for X'
util_ox = function(x, y, check = FALSE) {
  if (x$rk > 0) {
    out <- x$u[,seq_len(x$rk)]
    if (!missing(y))
      out <- out %*% y
    if (check)
      out[abs(out) < x$tol] <- 0
  } else
    out <- matrix(0, nrow(x$X), 1)
  return(out)
}

'orthonormal basis sets for t(X)'
util_oxp = function(x, y, check = FALSE) {
  if (x$rk > 0)
    x$v[, seq_len(x$rk)]
  else
    matrix(0, ncol(x$X), 1)
}

'pseudo inverse of X'
util_pinvx = function(x, y, check = FALSE) {
  if (x$rk > 0 ) {
    out <- x$v(, seq_len(x$rk)) %*% 
      tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]),
                 x$u[, seq_len(x$rk)])
  } else
    out <- matrix(0, ncol(x$X), nrow(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'pseudo-inverse of t(X)'
util_pinvXt = function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- x$u[, seq_len(x$rk)] %*% 
      tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]),
                 x$v[, seq_len(x$rk)])
  else
    out <- matrix(0, nrow(x$X), ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

'pseudo-inverse of crossprod(X)'
util_pinvXtX = function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- x$v[, seq_len(x$rk)] %*%
      tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

'pseudo-inverse of tcrossprod(X)'
util_-invXXt <- function(x, y, check = FALSE) {
  if (x$rk > 0) {
    out <- x$u[, seq_len(x$rk)] %*%
      tcrossprod(diag(x$d[seq_len(x$rk)]^(-2)),
                 x$u[, seq_len(x$rk)])
  } else
    out <- rep(0, nrow(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'computes u * (diag(x^n)) * t(u)'
util_power <- function(x, n, y, check = FALSE) {
  if (missing(n))
    n <- 1
  if (x$rk > 0) {
    if (missing(y))
      out <- x$u[, seq_len(x$rk)] %*% tcrossprod(diag(x$d[seq_len(x$rk)] ^ n), x$u[, seq_len(x$rk)])
    else
      out <- x$u[, seq_len(x$fk)] %*% diag(x$d[seq_len(x$fk)] ^ n) %*% crossprod(x$u[, seq_len(x$fk)], y)
  } else {
    if (missing(y))
      out <- rep(0, nrow(x$X))
    else
      matrix(0, nrow(x$X), ncol(y))
  }
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'computes v * (diag(x^n)) * t(v)'
util_powert <- function(x, n, y, check = FALSE) {
  if (missing(n))
    n <- 1
  if (x$rk > 0)
    out <- x$v[, seq_len(x$rk)] %*% tcrossprod(diag(x$d[seq_len(x$rk)]^(n)), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  
  return(out)
}

'residual projecting matrix'
util_r <- function(x) {
  diag(ncol(x$X)) - .self$getOP(x)
}

'computation of crossprod(X)'
util_xpx = function(x, y, check = FALSE) {
  if (x$rk > 0) {
    out <- x$v[, seq_len(x$rk)] %*%
      tcrossprod(diag(x$d[seq_len(x$rk)]^2),
                 x$v[, seq_len(x$rk)])
  } else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'computation of tcrossprod(X)'
util_xxp <- function(x, y, check = FALSE) {
  if (x$rk > 0) {
    if(missing(y)) {
      out <- x$u[, seq_len(x$rk)] %*%
        tcrossprod(diag(x$d[seq_len(x$rk)]^2),
                   x$u[, seq_len(x$rk)])
    } else
      x$u[, seq_len(x$rk)] %*% diag(x$d^2) %*%
      crossprod(x$u[, seq_len(x$rk)], y)
  } else
    out <-rep(0, nrow(x$X))
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

'set contrast fields'
set_con <- function(x, name, fieldType, c, iX0 = FALSE) {
  if (missing(x))
    x <- KWX
  out <- list()
  out$c <- matrix(c, ncol = 1)
  out$name <- name
  out$fieldType <- fieldType
  if (!iX0) {
    out$iX0 <- "c"
    if (length(c) == 0) {
      tmp <- .self$C2Tsp(x, c(), plus = TRUE)
    } else {
      tmp <- .self$C2Tsp(x, c, plus = TRUE)
    }
    out$X1$uKX1 <- tmp[[1]]
    out$X0$uKX0 <- tmp[[2]]
  } else {
    iX0 <- c
    iX0 <- .self$getiX0check(iX0, sL)
    out$iX0 <- iX0
    out$X0$uKX0 <- tcrossprod(.self$getOX(x), .self$getXi(x, iX0))
    if (length(iX0) == 0) {
      out$c <- .self$getXtX(x)
      out$X10$uKX1 <- .self$getCUKX(x)
    } else {
      out$c <- .self$geti02C(x, iX0)
      out$X1$uKX1 <- .self$getC2Tsp(x, out$c, plus = TRUE)
    }
  }
  
  C <<- lappend(C, out)
}

'traces for effective df calculation. If oneout = TRUE then returns trRV.'
util_trRV <- function(oneout = FALSE) {
  rk <- KWX$rk
  sL <- nrow(KWX$X)
  
  u <- KWX$u[, seq_len(rk)]
  if (oneout) {
    if (rk == 0)
      return(0)
    else
      trmv <- sum(u %*% (XV %*% u))
    return(sum(diag(XV)) - trmv)
  } else {
    if (rk == 0) {
      trmv <- 0
      tmp <- norm(XV, "f")^2
      trv <- sum(diag(XV))
    } else {
      Vu <- XV %*% u
      trv <- sum(diag(XV))
      tmp <- norm(XV, "F")^2
      tmp <- tmp - 2 * norm(Vu, "f")^2
      dims$trRVRV <<- tmp + norm(crossprod(u, Vu,), "f")^2
      trmv <- sum(u %*% Vu)
    }
    dims$trRV <<- trV - trmv
  }
  dims$rdf <<- (trRV^2) / trRVRV
  dims$npred <<- ncol(X)
}

'compute the traceof MV, MVMV, and find the degrees of interest'
util_trMV <- function(X1, xv, oneout = FALSE) {
  u <- X1$u[, seq_len(KWX$rk)]
  if (oneout) {
    out <- sum(crossprod(u, crossprod(u, xv)))
  } else {
    Vu <- xv %*% u
    out <- list()
    out$trMV <<- sum(u %*% Vu)
    out$trMVMV <<- norm(crossprod(u, Vu), "f")^2
    out$idf <<- (trMV^2) / trMVMV
  }
  return(out)
}

'return residual forming matrix or set residuals'
util_res <- function(check = FALSE) {
  KWY <- rftFilter(K, W %*% modData$imageMatrix)
  if (KWX$rk < (nrow(KWX$X) - KWX$rk)) {
    res <<- KWY - KWX$u[, seq_len(KWX$rk)] %*%
      crossprod(KWX$u[, seq_len(KWX$rk)], KWY)
  } else {
    if (ncol(y) < (5 * nrow(KWX$X)))
      res <<- .self$getR() %*% KWY
    else {
      n <- .self$getNp(KWX)
      res <<- n %*% crossprod(n, KWY)
    }
  }
  mrss <<- colSums(res^2) / dims$trRV
  if (control$scaleResid)
    resid <<- t(t(tmp) / as.numeric(mrss))
  if (check)
    res[abs(res < KWX$tol)] <<- 0
}

'set the filtered and whitened design matrix'
util_set <- function(x) {
  out <- list()
  out$X <- x
  svdx <- svd(t(x))
  out$v <<- svdx$v
  out$u <<- svdx$u
  out$d <<- svdx$d
  out$tol <<- max(dim(X)) * max(absx$d) * .Machine$double.eps
  out$rk <<- sum(x$d > x$tol)
  return(out)
  }

## from spm_SpUtil
util_iscon = function() {
  'test whether weight vectors specify contrast'
}

'space tested while keeping size of X$i0'
util_i02X1 = function(plus = FALSE) {
  c <- util_i02C(x, i0)
  if (plus)
    out <- c2tsp(x, c, plus = TRUE)
  else
    out <- c2tsp(x, c)
  return(out)
}

'get the estimable parts of C0 and C1'
util_i02C = function(x, i0) {
  sL <- ncol(x$X)
  i0 <- util_check_i0(i0, sL)
  
  c0 <- diag(sL)
  c0 <- c0[, i0)]
  c1 <- diag(sL)
  c1 <- c1[, setdiff(seq_len(sL), i0)]
  
  if (!util_isinspp(x, c0))
    c0 <- util_opp(x, c0)
  if (!util_isinspp(x, c1))
    c1 <- util_opp(x, c1)
  
  if (!(length(c1) == 0)) {
    if (!(length(c0 == 0)))
      out <- util_r(util_set(c0), c1, check = TRUE)
    else
      out <- util_XpX(x)
  } else
    out <- c()
  return(out)
}

util_X02C = function(x, y, plus = FALSE) {
  if (plus) {
    if (!is.list(cukX0))
      X0 <- c()
    else
      X0 <- util_OX(x) %*% cukX0
  }
  out <- util_X02C(X0, x)
  return(out)
}

'effective F degrees of freedom dof(idf, rdf)'
util_i02edf = function(x, i0, V) {
  if (util_isspc(x))
    x <- util_set(x)
  i0 <- util_check_i0(i0, ncol(x$X))
  if (misisng(V))
    V <- diag(nrow(x$X))
  
  r <- util_trMV(x, V)
  m <- util_trMV(util_i02X1(x, i0), V)
  
  return(c(m$trMV^2 / m$trMVMV, r$trMV^2 / r$trMVMV))
}

util_size = function() {
  'size of design matrix'
}

util_iX0check = function(i0, sL) {
  'check iX0'
  util_check_i0(i0, sL)
}

util_check_i0 = function(i0, sL) {
  'check i0'
  if (any(i0 == 1 | 0) && (length(i0) == sL))
    i0c <- (seq_len(sL)[i0 != 0])
  else if (all(dim(i0) > 0) && any(floor(i0) != i0) || any(i0 < 1) || any(i0 > sL))
    stop('logical mask or vector of column indices required')
  else
    i0c= i0
  return(i0c)
}

util_X0_2_c <- function(X0, x) {
  'get the orthogonal compliment and project onto the estimable space'
  if (is.null(X0)) {
    sc0 <- utilSet(utilX(x, X0))
    if (sc0$rk)
      c <- utiloPp(x, utilR(sc0))
    else
      c <- utiloPp(x)
    
    # dont know if this is equivalent to matlab command c  = c(:,any(c));
    c <- c[, colSums(abs(c)) > 0]
    sL <- ncol(x$X)
    
    if ((ncol(c) == 0) && (ncol(X0) != sL))
      c <- rep(0, sL)
  } else
    c <- util_xpx(x)
  return(c)
}

## from sp_FcUtil
## fcortho
util_X1 = function(Fc, x) {
  if (!is.null(Fc$X0))
    util_ox(x) %*% Fc$X1o$ukX1o
  else
    Fc$X1
}

util_X0 <- function(Fc, x) {
  if (!is.null(Fc$X0))
    util_ox(x) %*% Fc$X0$uKX0
  else
    Fc$X0
}

util_isnull = function(Fc, x) {
  any(colSums(abs(util_opp(x, Fc$c))) > 0)
}

util_is_T = function(x, c) {
  boul <- 1
  if (!util_isinspp, x, c)
    c <- util_opp(x, c, check)
  if ((rank(c) > 1) || any(crossprod(c) < 0))
    boul <- 0
  return(boul)
}

util_FconFields = function() {
  c("name", "fieldType", "c", "X0", "iX0", "X1", "idf", "Vcon", "Vspm")
}

util_minFc = function() {
  minFc <- list()
  minFC$name <- c()
  minFC$fieldType <- c()
  minFC$c <- c()
  minFC$X0 <- c()
  minFC$X1o <- c()
}

util_IsFcon = function(Fc) {
  b <- 1
  minnames <- fieldnames(sf_MinFcFields);
  FCnames <- fieldnames(Fc);
  for (i in 1:length(minnames)) {
    b <- (b && any(minnames[i] == FCnames))
    if !b
      break
  }
}

util_IsSet = function(Fc) {
  !util_isempty_X0(Fc) | !util_isempty_X1o(Fc)
}

util_isemptyX1 = function(Fc) {
  if (!is.null(Fc$X0)) {
    b <- is.null(Fc$X1o$ukX1o);
    if b != is.null(Fc$c), 
      Fc$c, Fc$X1o$ukX1o, error('Contrast internally not consistent')
  } else {
    b = is.null(Fc$X1o);
    if b != is.null(Fc$c), 
      Fc$c, Fc$X1o, error('Contrast internally not consistent');
  }
}

util_X1 = function(Fc, sX) {
  if (!is.null(Fc$X0))
    util_ox(sX) %*% Fc$X1o$uKX1
  else
    Fc$X1o;
}

util_X0 = function(Fc, sX) {
  if (!is.null(Fc$X0))
    util_ox(sX) %*% Fc$X0$uKX0; 
  else
    Fc$X0;
}

util_isempty_X0 = function(Fc) {
  if (!is.null(Fc$X0))    
    is.null(Fc$X0$ukX0); 
  else
    is.null(Fc$X0)
}

util_hsqr = function(Fc, sX) {
  if (!is.null(Fc$X0))
    crossprod(util_ox(util_set(Fc$X1$ukX1)), util_cukx(sX))
  else
    crossprod(util_ox(util_set(Fc$X1)), x$X)
}

util_H = function(Fc, sX) {
  if (!is.null(Fc$X0))
    crossprod(util_hsqr(Fc, sX))
  else
    Fc$c %*% tcrossprod(MASS::ginv(crossprod(Fc$X1)), Fc$c)
}

util_Y0 = function(Fc, x, b) {
  x$X %*% (diag(ncol(x$X)) - util_xpx(sX)) %*% util_H(Fc, x)) %*% b;
}

util_fcortho = function(Fc1, sX, Fc2) {
  c1_2  = Fc1$c - util_H(Fc2,sX) %*% util_xpx(sX,Fc1$c)
  
  Fc1o  = spm_FcUtil('Set', ['(' Fc1$name ' |_ (' Fc2.name '))'], Fc1$fieldType, 'c+', c1_2, sX);
}

util_Rortho = function (Fc1, sX, Fc2) {
  if (length(Fc2) == 0) {
    if (length(Fc1) <= 1)
      b <- 0
    else {
      c1 <- Fc1[[1]]$c
      for (i in 2:length(Fc1))
        c1 <- cbind(c1, Fc1[[i]]$c)
      b <- !any(abs(upper.tri(crossprod(c1, util_xpx(sX, c1)))) > sX$tol)
    }
  } else {
    c1 <- Fc1[[1]]$c
    for (i in 2:length(Fc1))
      c1 <- cbind(c1, Fc1[[i]]$c)
    c2 <- Fc2[[1]]$c
    for (i in 2:length(Fc2))
      c2 <- cbind(c1, Fc1[[i]]$c)
    b <- !any(abs(crossprod(c1, util_pinvxpx(x, c2)) > sX$tol))
  }
  return(b)
}

util_in = function(Fc1, x, Fc2) {
  l1 <- length(Fc1)
  l2 <- length(Fc2)
  for (j in seq_len(l1)) {
    if (!util_isinspp(x, Fc1[[j]]$c))
      c1 <- util_opp(x, Fc1[[j]]$c)
    else
      c1 <- Fc1[[j]]$c
    sc1 <- util_set(c1)
    S <- Fc1[[j]]$fieldType
    
    boul = 0
    idxFc1 <- c()
    idxFc2 <- c()
    for (i in seq_len(l2)) {
      if (Fc2[[i]]$fieldType == S) {
        boul <- util_equal(sc1, util_opp(x, Fc2[[i]]$c))
        
        if (boul && S == "T") {
          atmp <- util_x1o(Fc1[[j]], x)
          btmp <- util_x1o(Fc2[[i]], x)
          boul <- !any(crossprod(atmp, btmp) < 0)
        }
        if (boul) {
          idxFc1 <- c(idxFc1, j)
          idxFc2 <- c(idxFc2, i)
        }
      }
    }
  }
  return(list(idxFc1, idxFc2))
}

util_notunique = function(Fc, x) {
  l <- length(Fc)
  if (l == 1)
    out <- c()
  else {
    out <- list(1 + util_in(Fc[1], x, Fc[2:l]), 1 + util_notnunique(Fc[2:l], x))
    for (i in length(out))
      out[[i]] <- out[[i]] + 1
  }
  return(out)
}

setContrast = function(object, c, fieldType = "T", action = "c", plus = FALSE) {

}
