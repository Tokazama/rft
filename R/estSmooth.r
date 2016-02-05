#' Estimates smoothness of an image
#' 
#' The partial derivatives of an image in x, y, and z directions are used
#' to create a covariance matrix which in turn is used to calculate the full-widths
#' at half maxima (FWHM). The FWHM is equivalent to the estimated image smoothness.
#' 
#' The resels per voxel image (\code{RPVImg}) represents the estimated resel at each 
#' individual voxel throughout the search region. This is used in place of volumetric
#' measurements (or sum voxel measurements) when estimating the p-value of a cluster
#' using \code{rft.pval}. Using resels per voxel to estimate cluster level statistics
#' attempts to offset the natural probability to obtain significant clusters solely by
#" chance in very smooth regions at low thresholds.
#'
#' It is possible to use a single computed statistical parametric map image to
#' estimate the FWHM. However, it's recommended that FWHM estimates are obtained
#' from the residuals of statistical models (Stefan J.K et al., 1999). Therefore, this function
#' is optimized to estimate the pooled smoothness of the residual images from a linear model. 
#' 
#' A scaling factor is used to correct for differences when using the \code{sample} option.
#' Scaling isn't effective when the number of images is very low and typically results in an overestimation of the the FWHM. If only one image or numeric
#' vector is entered then the scaling factor isn't used. If a numeric vector is entered the \code{imageMake}
#' function is used to prepare it for smoothness estimation (see Worsley et al., 1999).
#'
#' @param img object of class antsImage
#' @param mask input mask, must match matrix
#' @param df degrees of freedom
#' @param makeRPV make resels per voxel image?
#' @param sample number of images to use for estimating smoothing
#' @return Outputs the estimated covariance matrix and fwhm in as an object of class list.
#' @references
#' Hayasaka (2004) Nonstationary cluster-size inference with random field and permutation methods.
#' Worsley K.J. (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worsley K.J. (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Worsley K.J. (1999) Detecting Changes in Nonisotropic Images
#' Stefan J.K. (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
#' @examples
#' # estimatation of a single images smoothness
#' outimg1 <-makeImage(c(10,10,10), rnorm(1000))
#' maskimg <-getMask(outimg1)
#' fwhm <-est.Smooth(outimg1,maskimg)
#' # estimation of smoothness of overall sample images in a statistical model
#' outimg2 <-makeImage(c(10,10,10), rnorm(1000))
#' imat <-imageListToMatrix(list(outimg1, outimg2), maskimg)
#' variable <-rnorm(2)
#' fit <-lm(imat ~ variable)
#' res <-residuals(fit)
#' 
#' @export estSmooth
estSmooth <-function(x,mask,df,makeRPV=TRUE,sample){
  D <-mask@dimension
  dimx <-1:dim(mask)[1]
  dimy <-1:dim(mask)[2]
  dimz <-1:dim(mask)[3]
  dimx1 <-2:(dim(mask)[1]+1)
  dimy1 <-2:(dim(mask)[2]+1)
  dimz1 <-2:(dim(mask)[3]+1)
  dimx2 <-3:(dim(mask)[1]+2)
  dimy2 <-3:(dim(mask)[2]+2)
  dimz2 <-3:(dim(mask)[3]+2)
  
  if (class(x)=="antsImage"){
    imgar <-as.array(x)
    classval <-1
    scale <-1
    n <-1
  }else if(class(x)=="numeric"){
    classval <-2
    x <-matrix(x,nrow=1)
    imgar <-as.array(makeImage(mask,x))
    scale <-1
    n <-1
  }else if(class(x)=="matrix"){
    classval <-3
    if (missing(sample)){
      nfull <-nrow(x) # original number of images (rows)
    }else{
      nfull <-nrow(x)
      rsamples <-sample(nrow(x),sample)
      x <-x[rsamples,]
    }
    n <-nrow(x) # number of images in sample (rows)
    scale <-(nfull/df)*(1/n)
  }
  # used to find partial derivatives (taken out of loop to save time)
  m<-d<-array(0, dim=dim(mask)+2)
  Vxx<-Vyy<-Vzz<-Vxy<-Vxz<-Vyz<-array(0,dim=dim(mask))
  
  maskar <-as.array(mask)
  m[dimx1,dimy1,dimz1] <-maskar
  xm <-m[dimx,dimy1,dimz1]+m[dimx1,dimy1,dimz1]+m[dimx2,dimy1,dimz1]
  ym <-m[dimx1,dimy,dimz1]+m[dimx1,dimy1,dimz1]+m[dimx1,dimy2,dimz1]
  zm <-m[dimx1,dimy1,dimz]+m[dimx1,dimy1,dimz1]+m[dimx1,dimy1,dimz2]
  
  progress <-txtProgressBar(min=0, max=n, style=3)
  for (i in 1:n){
    if (classval > 1){
      imgar <-as.array(makeImage(mask,x[i,]))
    }
    #calculate partial derivatives x
    d[dimx1,dimy1,dimz1] <-imgar
    dx <-(imgar-d[dimx,dimy1,dimz1])
    #calculate partial derivatives y
    dy <-(imgar-d[dimx1,dimy,dimz1])
    #calculate partial derivatives z
    dz <-(imgar-d[dimx1,dimy1,dimz])
    
    Vxx <-Vxx+(dx*dx)
    Vyy <-Vyy+(dy*dy)
    Vzz <-Vzz+(dz*dz)
    if (makeRPV=="TRUE"){
      Vxy <-Vxy+(dx*dy)
      Vxz <-Vxz+(dx*dz)
      Vyz <-Vyz+(dy*dz)
    }
    setTxtProgressBar(progress,i)
  }
  close(progress)
  rm(dx,dy,dz)
  Vxx <-Vxx*scale
  Vyy <-Vyy*scale
  Vzz <-Vzz*scale
  if (makeRPV=="TRUE"){
    Vxy <-Vxy*scale
    Vxz <-Vxz*scale
    Vyz <-Vyz*scale
    # calculated similarly to SPM. I believe it's derived from Worsley et al., (1999):
    # dif_u =(u1-u0,...,uD-u0)
    # Resels=(1/factorial(D))*sum(|dif_u'dif_u|^(1/2)*(4*log(2))^(-D/2))
    rpv <-Vxx*Vyy*Vzz+Vxy*Vyz*Vxz*2-Vyz*Vyz*Vxx-Vxy*Vxy*Vzz-Vxz*Vxz*Vyy
    rpv[rpv < 0] <-0
    rpv <-sqrt(rpv/(4*log(2))^D)
    # resels per voxel image
    rpv <-rpv*maskar # ensure no voxels exist outside mask
    RPVImg <-makeImage(mask,rpv)
  }
  xyzm <-xm+ym+zm
  xyzm[xyzm !=9] <-0
  xyzm[xyzm==9] <-1
  nvox <-sum(xyzm)
  xyz <-cbind((Vxx*xyzm),(Vyy*xyzm),(Vzz*xyzm))
  xyz <-sqrt(xyz/(4*log(2)))
  xyz <-colSums(xyz)/nvox
  if (makeRPV=="TRUE"){
    rpv <-sum(rpv)/nvox
    R <-rpv^(1/D)*(xyz/prod(xyz)^(1/D))
    fwhm <-1/R
    results <-list(fwhm=fwhm,RPVImg=RPVImg,rpv=rpv,xyz=xyz)
  }else if (classval==1 | classval==2){
    fwhm <-1/xyz
    results <-list(fwhm=fwhm,xyx=xyz)
  }
  return(results)
}
