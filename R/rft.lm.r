#' Fits a linear model according to random field theory
#' 
#' @param formula
#' @param conmat contrast matrix of contrasts by conditions (see example for more details)
#' @param mask object of type \code{antsImage}. Typically the same mask used to produce the image matrix.
#' @param tol (default = 1e-07)
#' @param statdir directory that statistical images are saved to (if not specified results are not saved)
#' @param sample how many of the residual images will be used to estimate FWHM and resels (default is all images)
#' @param findSmooth if \cod{findSmooth="FALSE"} then FWHM and resels aren't calculated (default is \code{"TRUE"}
#' @param concise if \code{concise="FALSE"} then residuals and coefficients are returned (default is \code{"TRUE"})
#' @param Intercept include intercept in design model
#' @param subMean scales each subject by image mean intensity (also uses a rough scale to physiologically meaningful values in PET (50 mls of blood / 100 mls of brain tissue /min))
#' @param regMean creates a control variable for each subjects mean image intensity value
#' @param groupMean creates a control variable for the mean image intensity value of each level of the identified factor 
#' @param grandMean scales the entire by the mean of all images
#' @param GG perform Greenhouse-Geisser correction for degrees of freedom?
#' @param verbose
#'
#' @return A list comprising:
#' \itemize{
#' \item{conmat}{contrast matrix used to create \code{Timgs}}
#' \item{Designmatrix}{the design matrix used to fit the model}
#' \item{df}{degrees of freedom}
#' \item{residuals}{image matrix of residuals}
#' \item{coefficients}{image matrix of coefficients}
#' \item{Timgs}{list of the t-field images } 
#' \item{fwhm}{the estimated FWHM}
#' \item{RPVImg}{resels per voxel image}
#' \item{resels}{estimated resels for the search space}
#' }
#' @description
#'
#' The images are each mean scaled to zero prior to
#' fitting (Friston K.J., 1995). As recommended by Stefan et al. (1999), the residuals are used to 
#' estimate the full width at half-maxima (FWHM; or "smoothness") of the final statistical image. The
#' residuals are standardized using the mean sum of squares of the residuals.
#' 
#' Standardized Residuals = residuals/sqrt(MRSS/Error of Degrees of Freedom)
#' 
#' The the resolutions in voxel space (resels) are calculated using the estimated FWHM and the image search 
#' region (mask). The StatImgs output contains a list of t-statistical fields in image space
#' determined according to specified contrasts (\code{conmat}).
#' 
#'
#'
#' @References 
#' Friston K.J., (1995) Statistical Parametric Maps in Functional Imaging: A General Linear Approach 
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worlsey K.J., et al. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Stefan J.K., (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' 
#' @Author Zachary P. Christensen
#' @kewords rft.pval, rft.euler, rft.resels
#' @note: function currently in beta phase. Waiting for acceptance of peer-reviewed paper
#' @examples
#' 
#'
#' @export rft.lm
rft.lm <-function(formula, conmat, conType, mask, 
                  tol=1e-07, statdir=NULL, sample=NULL, 
                  intercept=TRUE, subMean=TRUE, regMean=FALSE, 
                  groupMean=NULL, findSmooth=TRUE, concise=TRUE, 
                  grandMean=FALSE, verbose=TRUE){
  
  Y <-model.frame(formula)[[1L]]
  nsub <-nrow(Y)
  nvox <-sum(as.array(mask))
  #myresp <-as.character(formula(formula)[2])
  
  if (subMean=="TRUE" | regMean=="TRUE" | !is.null(groupMean)){
    g <-rowSums(Y)/nvox
    
    if (subMean=="TRUE"){
      X <-model.matrix(update(formula, ~ . - 1))
      for (sub in 1:nsub){
        Y[nsub,] <-Y[nsub,]*(50/g[nsub])
        X[nsub,] <-X[nsub,]*(g[nsub]/50)
      }
    }
    if(regMean=="TRUE"){
      X <-model.matrix(update(formula, ~ . + g - 1))
      conmat <-cbind(conmat,0)
    }
    if(!is.null(groupMean)){
      gg <-ave(g, groupMean)
      X <-model.matrix(update(formula, ~ . + gg - 1))
      conmat <-cbind(conmat,0)
    }
  }else{
    X <-model.matrix(update(formula, ~ . - 1))
  }
  # Grand Mean Scaling
  if (grandMean=="TRUE")
    Y <-Y * mean(Y[Y !=0])
  
  if (intercept=="TRUE"){
    X <-cbind(X, 1L)
    colnames(X)[ncol(X)] <-"Intercept"
    conmat <-cbind(conmat,0)
    colnames(conmat)[ncol(conmat)] <-"Intercept"
  }
  
  #### Fit General Linear Model ####
  # fit using OLS and calculate contrasts
  z <-rft.fit(X, Y, conmat, conType, statdir=statdir)
  
  #### Estimate FWHM/RESELS ####
  if (findSmooth=="TRUE"){
    if (verbose)
      cat("Estimating FWHM. \n")
    
    sr <-z$residuals/sqrt(z$rss/z$df[2])
    fwhm <-estSmooth(sr, mask, df, sample=sample)
    ans <-lappend(ans,fwhm$fwhm)
    names(ans)[length(ans)] <-"fwhm"
    ans <-lappend(ans,fwhm$RPVImg)
    names(ans)[length(ans)] <-"RPVImg"
    
    if (verbose)
      cat("Calculating resels. \n")
    
    resels <-rft.resels(mask, fwhm$fwhm)
    ans <-lappend(ans,resels)
    names(ans)[length(ans)] <-"resels"
  }
  
  #### Prepare Output ####
  ans <-list(X, resels, fwhm$fwhm, fwhm$RPVImg, z$coefficients, z$qr, z$df, z$StatImgs)
  names(ans) <-c("design","resels","fwhm","RPVImg","coefficients","qr","df","StatImgs")
  
  if (concise=="FALSE"){
    ans <-lappend(ans, z$residuals)
    names(ans)[length(ans)] <-"residuals"
  }
  ans
}
