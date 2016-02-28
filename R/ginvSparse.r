#' Generalized Inverse of a Matrix (sparse)
#'
#' Identical to the \code{ginv} function from the MASS package but uses \code{svdSparse} instead of \code{svd}.
#'
#' @param X Matrix for which the Moore-Penrose inverse is required
#' @param tol relative tolerance to detect zero singular values
#'
#' @seealso \code{ginv} , \code{svd}, \code{svdSparse}
#' @examples
#' 
#' X <-Matrix::Matrix(rnorm(100),10,10,sparse=TRUE)
#' Xinv <-ginvSparse(X)
#'
#' @export ginvSparse
ginvSparse <-function(X, tol=.Machine$double.eps){
    Xsvd <- svdSparse(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}
