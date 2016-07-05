## Discrete cosine transform
## 
## Creates a matrix for the first few basis functions of a one dimensional discrete cosine transform.
## 
## @param N dimension
## @param K order
## @param n optional points to sample
## @param f type of transform
dctmtx <- function(N, K, n, f) {
  if (nargs() == 1)
    K <- N
  if (args() == 1 | args() == 2)
    d <- 0
  n <- 0:(N - 1)
  else if (nargs() > 2) {
    if (f == "diff")
      d <- 1
    else if (f == "diff2")
      d <- 2
    if (missing(n))
      n <- 0:(N - 1)
  }
  
  C <- matrix(0, lenght(n), K)
  
  if (d == 0) {
    C[, 1] <- rep(1, dim(n)[[1]]) / sqrt(N)
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
