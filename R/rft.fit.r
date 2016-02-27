#' Fits model uing RFT contrasts
#' 
#' 
#' 
#' @param X design matrix
#' @param Y image matrix
#' @param conmat matrix of contrasts
#' @param conType a list of specifying the type of statistical field to be computed (i.e. t-values="T" or F-values="F")
#'               corresponding to each row of conmat specifying of
#' @param statdir directory to save contrast images
#'
#'
#' @export rft.fit
rft.fit <-function(X, Y, conmat, conType, statdir){
  if (require("Matrix")){
    cat("Using sparse matrices for analysis. \n")
  }else {
    stop("Matrix package is required for this function. \n")
  }
  if (ncol(conmat) !=ncol(X)){
    stop("Object conmat must have the same number of columns as factors/variables included in the model.")
  }
  n <-nrow(X) # number of subjects

  z <-.lm.fit(X,Y,tol=tol)
  #  # TEST THIS: may not be worth the initial computational time investment
  #  z$qr <-Matrix(z$qr, sparse=TRUE)
  #  z$coefficients <-Matrix(z$coefficients, sparse=TRUE)
  #  z$residuals <-Matrix(z$residuals, sparse=TRUE)
  #  z$effects <-Matrix(z$effects, sparse=TRUE)
  #  
  #  
  
  
  rss <-colSums(z$residuals^2)
  Ip <-Diagonal(z$rank)
  In <-Diagonal(n)
  XX <-chol2inv(qr.R(x$qr))
  
  # df[degrees of interest, degrees of error]
  df <-c(z$rank-1,n-z$rank)

  # residual projecting matrix
  R <-In-X %*% ginvSparse(X)

  # retain value for computing F contrasts
  YRY <-t(Y) %*% R %*% Y

  qr <-z[c("qr", "qraux", "pivot", "tol", "rank")]
  # Q <-qr.Q(qr)
  
  StatImgs <-list()
  for (i in 1:nrow(conmat)){
    if (conType[i]=="T"){
      se <-t(as.matrix(sqrt((rss/df[2])*(conmat[i,] %*% XX %*% conmat[i,]))))
      se[se==0] <-.01 # control for NULLS that result from dividing by 0
      statmat <-(conmat[i,] %*% z$coefficients)/se
    }else if(conType[i]=="F"){
      c <-matrix(conmat[i,],ncol=1)
      X0 <-X %*% (Ip-c %*% ginvSparse(c))
      R0 <-In-X0 %*% ginvSparse(X0) # contrast projecting matrix
      M <-R0-R
      statmat <-((t(z$coefficients) %*% t(X) %*% M %*% z$coefficients)/z$YRY)*(df[1]/z$df[2])
    }
    statimg <-makeImage(mask, statmat)
    if (!missing(statdir)){
      antsImageWrite(statimg, file=paste(statdir, rownames(conmat)[i], ".nii.gz", sep=""))
    }
    StatImgs <-lappend(StatImgs,statimg)
    names(StatImgs)[length(StatImgs)] <-rownames(conmat)[i]
  }
  if (!is.null(rownames(conmat)))
    names(StatImgs) <-rownames(conmat)
  c(z[c("coefficients", "residuals", "effects")], list(qr = structure(qr, class = "qr"), df=df, StatImgs=StatImgs))
  }
