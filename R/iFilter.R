#' Filter 
#'
#' @param Y data matrix
#' @param K filter matrix
#' @param K struct array containing partition specific specifications
#' @param K[[s]]$RT observation interval in seconds
#' @param K[[s]]$row row of Y constituting block/partitions
#' @param K[[s]]$HParams cut-off period in seconds
#' @param K[[s]]$X0 low frequencies to be removed (DCT)
#' @return filtered data
#' 
#' out <- iModelMake(X = z$X, y = z$y[i], data = z$iData)
#'
#' @export iFilter
iFilter <- function(K, Y) {
  if (missing(Y) && is.data.frame(K)) {
    # set K$X0
    out <- list()
    for (s in seq_len(length(K))) {
      # create filter index from filter data.frame
      subout <- list()
      subout$row <- seq_len(nrow(K))[K$Filter == paste("F", s, sep = "")]
      subout$RT <- K$RT[subout$row[1]]
      subout$HParam <- K$HParam[subout$row[1]]
      
      # determine low frequencies to be removed
      nk <- length(subout$row)
      n <- floor(2 * (nk * subout$RT) / subout$HParam + 1)
      X0 <- .dctmtx(nk, n)
      subout$X0 <- X0[, 2:ncol(X0)]
      out[[i]] <- subout
    }
    return(out)
  } else {
    if (is.data.frame(K)) {
      # ensure requisite fields are present
      if (is.null(K[[1]]$X0))
        K <- iFilter(K)
      
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