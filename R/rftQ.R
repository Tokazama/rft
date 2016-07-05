## returns an (n x n) (inverse) autocorrelation matrix for an AR(p) process
##
## a  - vector of (p) AR coefficients
## n  - size of Q
## q  - switch to return inverse autocorrelation or precision [default q = 0]
## spm_Q uses a Yule-Walker device to compute K where:
## 
## y = K*z
## 
## such that y is an AR(p) process generated from an i.i.d innovation 
## z.  This means
## 
## cov(y) = <K*z*z'*K> = K*K'
## 
## If called with q ~= 0, a first order process is assumed when evaluating
## the precision (inverse covariance) matrix; i.e., a = a(1)
rftQ <- function(a, n, q) {
  # default
  try, q; catch, q = 0; end

if q

# compute P (precision)
  A    = [-a(1) (1 + a(1)^2) -a(1)];
Q    = spdiags(ones(n,1)*A,(-1:1),n,n);

else
  
  # compute Q (covariance)
  p <- length(a);
  A <- [1 -a(:)']
        P    = spdiags(ones(n,1)*A,-(0:p),n,n);
        K    = inv(P);
        K    = K.*(abs(K) > 1e-4);
        Q    = K*K';
        Q    = toeplitz(Q(:,1));
        
        end
}