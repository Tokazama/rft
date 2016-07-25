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