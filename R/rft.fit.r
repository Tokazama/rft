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
  rss <-colSums(z$residuals^2)
  Ip <-Diagonal(z$rank)
  In <-Diagonal(n)
  XX <-chol2inv(qr.R(x$qr))
  
  # df[degrees of interest, degrees of error]
  z <-lappend(z,c(z$rank-1,n-z$rank))
  names(z)[length(z)] <-"df"
  
  # residual projecting matrix
  z <-lappend(z,In-X %*% pinv(X)) 
  names(z)[length(z)] <-"R"
  
  # retain value for computing F contrasts
  z <-lappend(z,t(Y) %*% R %*% Y) 
  names(z)[length(z)] <-"YRY"
  
  StatImgs <-list()
  for (i in 1:nrow(conmat)){
    if (conType[i]=="T"){
      se <-t(as.matrix(sqrt((rss/df[2])*(conmat[i,] %*% XX %*% conmat[i,]))))
      se[se==0] <-.01 # control for NULLS that result from dividing by 0
      statmat <-(conmat[i,] %*% z$coefficients)/se
    }else if(conType[i]=="F"){
      c <-matrix(conmat[i,],ncol=1)
      X0 <-X %*% (Ip-c %*% pinv(c))
      R0 <-In-X0 %*% pinv(X0) # contrast projecting matrix
      M <-R0-z$R
      statmat <-((t(z$coefficients) %*% t(X) %*% M %*% z$coefficients)/z$YRY)*(df[1]/z$df[2])
    }
    statimg <-makeImage(mask, statmat)
    if (!missing(statdir)){
      antsImageWrite(statimg, file=paste(statdir, rownames(conmat)[i], ".nii.gz", sep=""))
    }
    StatImgs <-lappend(StatImgs,statimg)
    names(StatImgs)[length(StatImgs)] <-rownames(conmat)[i]
  }
  z <-lappend(z,StatImgs)
  names(z)[length(z)] <-"StatImgs"
  
  z
  }
