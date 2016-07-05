## return error covariance constraints for basic ANOVA designs
##
## required fields:
## xVi.I    - n x 4 matrix of factor level indicators
##              I(n,i) is the level of factor i for observation n
## xVi.var  - 1 x 4 vector of flags
##              var(i) = 1; different variance among levels of factor i
## xVi.dep  - 1 x 4 vector of flags
##              dep(i) = 1;      dependencies within levels of factor i
##
## Output:
## xVi.Vi   -  cell of covariance components
## or
## xVi.V    -  speye(n,n)
##
nonSphericity <- function(xVi) {
  # create covariance components Q{:}
  n <- nrow(xVi$I)
  f <- ncol(xVi$I)
  l <- max(xVi$I)
      
  # if var(i): add variance component for each level of factor i,
  Q <- list()
  for (i in find(xVi.var)) {
    for (j in 1:l[i]) {
      u <- xVi$I[, i] == j
      q <- diag(u)
      Q <- lappend(Q, q)
    }
  }
  
  # effects (discounting factors with dependencies) as defined by interactions
  X <- rep(1, n)
  for (i in find(~xVi.dep && (l > 1))) {
    Xi <- sparse(1:n, xVi.I(:,i), 1, n, l[i])
    Xj <- X
    X <- sparse(n, 0)
    for (j in 1:ncol(Xi)) {
      for (k in 1:ncol(Xj)) {
        X[, end + 1] <- Xi[, j] & Xj[, k]
      }
    }
  }
        
  # dependencies among repeated measures created by the Hadamard product
  for (i in find(xVi$dep)) {
    q <- sparse(seq_len(n), xVi.I[, i], 1, n, l[i])
    P <- tcrossprod(q)
      for (j in 1:ncol(X)) {
        for (k in (j + 1):ncol(X)) {
          Q = lappend(Q, (tcrossprod(X[, j], X[, k]) +
                          tcrossprod(X[, k], X[, j])) %*% P)
        }
      }
  }
  
  # set Q in non-sphericity structure
  
  # if i.i.d nonsphericity (V) is known otherwise there are components {Vi}
  if (length(Q) > 1)
    xVi$Vi <- Q
  else
    xVi$V <- matrix(nrow = n, nrow = n)
}
