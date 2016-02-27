#' Singular Decomposition Using Sparse Matrices
#'
#'  
#'
#' @param x a matrix
#' @param nu number of left singular vectors to be computed.
#' @param nv number of right singular vectors to be computed.
#' 
#' @return
#' d a vector containing th esingular values of x
#' u a matrix whose columns contain the left singular ectors of x
#' v a matrix whose columns contain the right singular vectors of x
#' 
#' @details 
#' If x has less than three rows or columns then \code{La.svd} from the stats package is used. Otherwise \code{svds} from the aRPACK package is used.
#' 
#' @examples
#' 
#' 
#'
#' @export svdSparse
svdSparse <-function (x, nu=min(n,p), nv=min(n,p), k=min(n,p)){
  if (!usePkg("Matrix") | !usePkg("aRPACK")){
    stop("Need Matrix and aRPACK packages")
  }
  if (length(dim(x)) !=2) 
    stop("x must be a matrix")
  if (!is.matrix(X))
    X <- Matrix::Matrix(x,sparse=TRUE)
  dx <-dim(x)
  n <-dx[1L]
  p <-dx[2L]
  if (n >= 3 && p >= 3){
    z <-aRPACK::svds(x, k=k, nu=nu, nv=nv)
  }else{
    La.res <-La.svd(x, nu, nv)
    res <- list(d = La.res$d)
    if (nu) 
        res$u <- La.res$u
    if (nv) {
      if (is.complex(x)) 
        res$v <- Conj(t(La.res$vt))
      else res$v <- t(La.res$vt)
    }
    z <-list(Matrix(res$d,sparse=TRUE),Matrix(res$u,sparse=TRUE),Matrix(res$v,sparse=TRUE))
  }
  z
  }
