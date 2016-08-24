#' Restricted Maximum Likelihood for image class objects 
#' 
#' @param YY sample covariance matrix tcrossprod(Y) {Y = image matrix}
#' @param X design matrix
#' @param Q covariance components
#' @param N number of samples (default 1)
#' @param D flage for positive-definite scheme (default 0)
#' @param t regularisation (default 4)
#' @param hE hyperprior (default 0)
#' @param hP hyperprecision (default 1e-16)
#' @param its maximum iterations
#'
#' @return
#' \item{V} {m x m estimated errors = h[1]*Q[[1]] + h[2]*Q[[2]] + ...}
#' \item{h} {q x 1 ReML hyperparameters}
#' \item{Ph} {q x q conditional precision of h}
#' \item{F} {free energy F = log evidence = p(Y|X,Q) = ReML objective}
#' \item{Fa} {accuracy}
#' \item{Fc} {complexity (F = Fa - Fc)}
#' adapted from spm_ReML
#' @export iREML
iREML <- function(YY, X, Q, N, D, t, hE, hP, its) {
  if (missing(N))
    N <- 1
  if (missing(D))
    D <- 0
  if (missing(t))
    t <- 0
  if (missing(hE))
    hE <- 0
  if (missing(hP))
    hP <- 1e-16

  # ortho-normalise X----
  if (missing(X))
    stop("Please specify X")
  else
    X <- svd(X)$u
  if (!is.list(Q))
    Q <- list(Q)

  # dimensions----
  n <- nrow(Q[[1]])
  m <- length(Q)

  # catch NaNs----
  W <- Q
  #q <- is.finite(YY)
  #YY <- YY[q, q]
  #for (i in 1:m)
  #  Q[[i]] <- Q[[i]][q, q]

  # initialise h and specify hyperpriors----
  h <- matrix(0, m, 1)
  for (i in 1:m)
    h[i, 1] <- sum(diag(Q[[i]]))
  hE <- matrix(0, m, 1) + hE
  hP <- matrix(0, m, m) * hP
  dF <- Inf
  dh <- matrix(0, m, 1)
  dFdh <- matrix(0, m, 1)
  dFdhh <- matrix(0, m, m)
  PQ <- list()

  # ReML (EM/VB)----
  for (it in seq_len(its)) {
    # compute current estimate of covariance----
    #C <- matrix(0, n, n)
    C <- 0
    for (i in 1:m)
      C <- C + Q[[i]] * h[i]

    # positive [semi]-definite check----
    # might be able to use nlme::pdMat()
    for (i in 1:D) {
      if (min(eigen(C)$values) < 0) {
        t <- t - 1
        h <- h - dh
        dh <- .dx(dFdhh, dFdh, t)
        h <- h + dh
        C <- matrix(0, n, n)
        for (j in 1:m)
          C <- C + Q[[i]] * h[i]
      }
    }

    # E-step: conditional covariance cov(B|y) {Cq}----
    iC <- solve(C)
    iCX <- iC %*% X
    Cq <- solve(crossprod(X, iCX))
    # if (!any(X != 0))
    #   Cq <- solve(crossprod(X, iCX))
    # else
    #   Cq <- 0

    # M-step: ReML estimate of hyperparameters----
    P <- iC - iCX %*% Cq %*% t(iCX) # P = iV - iV %*% X %*% solve(crossprod(X, iV) %*% X) %*% crossprod(X, iV)
    U <- diag(n) - P %*% YY / N

    # dF/dh
    for (i in 1:m) {
      PQ[[i]] <- P %*% Q[[i]]
      dFdh[i] <- -sum(diag(PQ[[i]] %*% U)) * N / 2
    }

    # expected curvature E{dF / dhhh}
    for (i in 1:m) {
      for (j in 1:m) {
        # dF/dhh
        dFdhh[i, j] <- -sum(diag(PQ[[i]] %*% PQ[[j]])) * N / 2
        dFdhh[j, i] <- dFdhh[i, j]
      }
    }

    # add hyperpriors
    e <- h - hE
    dFdh <- dFdh - hP %*% e
    dFdhh <- dFdhh - hP

    # fisher scoring: update dh = -inv(ddF/dhh) * dF / dh
    dh <- .dx(dFdhh, dFdh, {t})
    h <- h + dh

    # predicted change in F - increase regularisation if increasing
    pF <- crossprod(dFdh, dh)
    if (pF > dF)
      t <- t - 1
    else
      t <- t + 1/4

    # final estimate of covarience (with missing data points)
    if (dF < 1e-1)
      break
  }

  # rebuild predicted covariance----
  V <- 0
  for (i in 1:m) {
    V <- V + W[[i]] * h[i]
  }

  # check V is positive semi-definite
  if (!D) {
    if (min(eigen(V)$values) < 0)
      iREML(YY, X, Q, N, 1, 2, hE[1], hP[1])
  }

  # log evidence = ln p(y|X, Q) = ReML objective = F = trace(R' *iC * R * YY) / 2 ...
  Ph <- - dFdhh


  if (nargs() > 4) {
    # tr(hP * inv(Ph)) - nh + tr...
    Ft <- sum(diag(hP %*% MASS::ginv(Ph))) - length(Ph) - length(Cq)

    # complexity
    Fc <- Ft / 2 +
      crossprod(e, hP) %*% e/2 +
      determinant(Ph %*% MASS::ginv(hP), logarithm = TRUE)$modulus / 2

    # accuracy - lnp(Y|h)
    Fa = Ft / 2 -
      sum(diag(C * P * YY * P)) / 2 -
      N * n * log(2 * pi) / 2 -
      N * determinant(C, logarithm = TRUE)$modulus / 2

    # free-energy
    FE <- Fa - Fc
    return(list(V = V, h = h, Ph = Ph, FE = FE, Fa = Fa, Fc = Fc))
  } else {
    return(list(V = V, h = h, Ph = Ph))
  }
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