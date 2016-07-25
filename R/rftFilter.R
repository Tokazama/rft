## @param Y data matrix
#'
#'
#' @param K filter matrix
#' @param K[[s]] struct array containing partition specific specifications
#' @param K[[s]]$RT observation interval in seconds
#' @param K[[s]]$row row of Y constituting block/partitions
#' @param K[[s]]$HParams cut-off period in seconds
#' @param K[[s]]$X0 low frequencies to be removed (DCT)
## 
#' @return filtered data
#'
rftFilter <- function(K, Y) {
  if (missing(Y) && is.list(K)) {
    # set K$X0
    for (s in seq_len(length(K))) {
      nk <- length(K[[s]]$row)
      n <- floor(2 * (nk * K[[s]]$RT) / K[[s]]$HParam + 1)
      X0 <- .dctmtx(nk, n)
      K[[s]]$X0 <- X0[, 2:ncol(X0)]
    }
    return(K)
  } else {
    if (is.list(K)) {
      # ensure requisite fields are present
      if (is.null(K[[1]]$X0))
        K <- rftFilter(K)
      
      # apply high pass filter
      if (length(K) == 1 && length(K[[1]]$row == ncol(Y)))
        Y <- Y - K[[1]]$X0 %*% crossprod(K[[1]]$X0, Y)
      else {
        for (i in seq_len(length(K))) {
          
          # select data
          y <- Y[K[[i]]$row, ]
          
          # apply high pass filter
          y <- y - K[[i]]$X0 %*% crossprod(K[[i]]$X0, y)
          
          # reset filtered data in Y
          Y[K[[i]]$row, ] <- y
        }
      }
    } else {
      Y <- K * Y
    }
    return(Y)
  }
}

# Discrete cosine transform
# 
# Creates a matrix for the first few basis functions of a one dimensional discrete cosine transform.
# 
# @param N dimension
# @param K order
# @param n optional points to sample
# @param f type of transform
.dctmtx <- function(N, K, n, f) {
  d <- 0
  
  if (nargs() == 1)
    K <- N
  if (args() == 1 | args() == 2)
    n <- seq_len(N - 1) - 1
  else if (nargs() > 2) {
    if (f == "diff")
      d <- 1
    else if (f == "diff2")
      d <- 2
    if (missing(n))
      n <- seq_len(N - 1) - 1
  }
  
  C <- matrix(0, lenght(n), K)
  
  if (d == 0) {
    C[, 1] <- rep(1, dim(n)[1]) / sqrt(N)
    for (k in 2:K)
      C[, k] <- sqrt(2 / N) * cos(pi * (2 * n + 1) * (k - 1) / (2 * N))
  } else if (d == 1) {
    for (k in 2:K)
      C[, k] = -2^(1 / 2) * (1 / N)^(1 / 2) * sin(1 / 2 * pi * (2 * n * k - 2 * n + k - 1) / N) * pi * (k - 1) / N;
  } else if (d == 2) {
    for (k in 2:K)
      C[, k] = -2^(1 / 2) * (1 / N)^(1 / 2) * cos(1 / 2 * pi * (2 * n + 1) * (k - 1) / N) * pi^2 * (k - 1)^2 / N^2;
  }
  return(C)
}

# @param dfdx = df/dx
# @param f = dx/dt
# @param t = integration time: (default t = Inf);
# if t is a cell (i.e., {t}) then t is set to:
# exp(t - log(diag(-dfdx))
#
.dx <- function(dfdx, f, t = Inf) {
  # t is a regulariser
  if (length(t) == 1)
    t <- exp(t - log(diag(-dfdx)))
  
  if (min(t) > exp(16)) {
    dx = - MASS::ginv(dfdx) %*% as.matrix(f, ncol = 1)
    dx =  as.array(dx, dim = dim(f));
  } else {
    # ensure t is a scalar or matrix
    if (length(t) > 1)
      t = diag(t)
    
    q <- matrix(0, nrow = max(dim(dfdx)) + 1, ncol = 1)
    q[1, 1] <- 1
    
    # augment Jacobian and take matrix exponential
    Jx <- rbind(0, cbind(t %*% f, t %*% dfdx))
    dx <- Matrix::expm(Jx) %*% q
    dx <- dx[2:nrow(dx), ]
  }
  dx
}
