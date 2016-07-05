## return error covariance constraints (for serially correlated data)
## @param v (1 x n) v(i) = number of observations for i-th block
## @param a AR coefficient expansion point  (default a = [])
## 
## a  = [] (default) - block diagonal identity matrices specified by v:
##
##   C{i}  = blkdiag( zeros(v(1),v(1)),...,AR(0),...,zeros(v(end),v(end)))
##   AR(0) = eye(v(i),v(i))
##
## otherwise:
##
##   C{i}     = AR(a) - a * dAR(a) / da;
##   C{i + 1} = AR(a) + a * dAR(a) / da;
##
## See also: rftQ
rftCe <- function(v, a) {
  if (missing(a))
    a <- c() 
  # create block diagonal components
    C <- list()
    l <- length(v)
    n <- sum(v)
    k <- 0
    if (l > 1) {
      for (i in seq_len(l)) {
        dCda <- rftCe(v[i], a)
        for (j in seq_len(length(dCda))) {
          [x,y,q] <- find(dCda{j});
           x <- x + k
           y <- y + k
           C <- lappend(C, sparse(x, y, q, n, n))
        }
        k <- v[i] + k
      }
    } else {
      # dCda
      if ~isempty(a) {
        Q    = spm_Q(a,v);
        dQda = spm_diff('spm_Q', a, v, 1);
        C[[1]] = Q - dQda[1] %*% a;
        C[[2]] = Q + dQda[1] %*% a;
      } else
        C[[1]] = speye(v,v);
    }
  return(C)
}
