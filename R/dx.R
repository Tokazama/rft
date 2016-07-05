## @param dfdx = df/dx
## @param f = dx/dt
## @param t = integration time: (default t = Inf);
## if t is a cell (i.e., {t}) then t is set to:
## exp(t - log(diag(-dfdx))
## 
## 
##
## 
dx <- function(dfdx, f, t = Inf) {
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
