#' Estimates smoothness of an image
#' 
#' The partial derivatives of an image in x, y, and z directions are used
#' to create a covariance matrix which in turn is used to calculate the full-widths 
#' at half maxima (FWHM). The FWHM is equivalent to the estimated image smoothness.
#'
#' It is possible to use a single computed statistical parametric map image to
#' estimate the FWHM. However, it's recommended that FWHM estimates are obtained
#' from the residuals of statistical models (Stefan J.K et al., 1999). 
#'
#' @param img - object of class antsImage
#' @param mask - input mask, must match matrix
#' @return Outputs the estimated covariance matrix and fwhm in as an object of class list.
#' @references
#' Hayasaka et al., (2004) Nonstationary cluster-size inference with random field and permutation methods.
#' Worsley K.J., (1992) A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain.
#' Worlsey K.J., (1996) A Unified Statistical Approach for Determining Significant Signals in Images of Cerebral Activation.
#' Stefan J.K., (1999) Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residual from the General Linear Model
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
#' DegFree <-fit$df.residual
#' meanres <-colMeans(res)
#' psd <-sqrt(colSums(res^2)/DegFree)
#' Zmat<-matrix(nrow=2,ncol=ncol(imat))
#' fwhm<-matrix(0L,nrow=1,ncol=3)
#' subs<-nrow(res)
#' for (i in 1:subs){
#' 	Zmat[i,]<-(res[i,]-Mres[1])/psd
#' 	img<-makeImage(mask,Zmat[i,])
#' 	smooth<-est.Smooth(img,mask)
#'	fwhm<-fwhm+smooth[[2]]
#' 	}
#' fwhm2<-sqrt(4*log(2)/(fwhm/(subs-1)))
#' 
#' @export estSmooth
estSmooth <-function(x,mask,df,sample){
  D <-mask@dimension
  dimx <-dim(mask)[1]
  dimy <-dim(mask)[2]
  dimz <-dim(mask)[3]
  if (class(x)=="antsImage"){
    classval <-1
    scale <-1
    n <-1
  }else if(class(x)=="numeric"){
    classval <-2
    x <-matrix(x,nrow=1)
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
    scale <-(nfull/(df[2]))*(1/n)
  }
  maskar <-as.array(mask)
  
  # used to find partial derivatives (taken out of loop to save time)
  d1<-dx2<-dy2<-dz2<-m1<-xm2<-ym2<-zm2<-Vxx<-Vyy<-Vzz<-Vxy<-Vxz<-Vyz<-array(0, dim=dim(mask)+1)
  maskar <-as.array(mask)
  m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
  xm2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
  ym2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
  zm2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
  xm3 <-(m1 + xm2)==2
  ym3 <-(m1 + ym2)==2
  zm3 <-(m1 + zm2)==2
  m3 <-(xm3*ym3*zm3)
  
  progress <-txtProgressBar(min=0, max=n, style=3)
  for (i in 1:n){
    if (classval > 1){
      img <-makeImage(mask,x[i,])
      }
    imgar <-as.array(img)
    #calculate partial derivatives x
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
    dx2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
    dx <-(d1-dx2)*m3
    #calculate partial derivatives y
    dy2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
    dy <-(d1-dy2)*m3
    #calculate partial derivatives z
    dz2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
    dz <-(d1-dz2)*m3
    Vxx <-Vxx+(dx*dx)
    Vyy <-Vyy+(dy*dy)
    Vzz <-Vzz+(dz*dz)
    Vxy <-Vxy+(dx*dy)
    Vxz <-Vxz+(dx*dz)
    Vyz <-Vyz+(dy*dz)

    setTxtProgressBar(progress,i)
  }
  close(progress)
  rm(dx,dy,dz)
  
  dx <-Vxx*scale
  dy <-Vyy*scale
  dz <-Vzz*scale
  dxy <-Vxy*scale
  dxz <-Vxz*scale
  dyz <-Vyz*scale
  
  xyz <-cbind(dx*m3,dy*m3,dz*m3)
  resel.img <-(dx*dy*dz)+(dxy*dyz*dxz*2)-(dx*dyz*dyz)-(dxy*dxy*dz)-(dxz*dy*dxz)  
  resel.img[resel.img < 0] <-0
  resel.img <-sqrt(resel.img/(4*log(2))^D)
  xyz <-sqrt((xyz)/(4*log(2)))
  
  # resels per voxel
  RPV <-resel.img[2:(dimx+1),2:(dimy+1),2:(dimz+1)]
  RPV <-RPV*maskar
  RPVImg <-as.antsImage(RPV)
  
  nvox <-sum(m3)
  resel.img <-sum(resel.img)/nvox
  xyz <-colSums(xyz)/nvox
  resels <-resel.img^(1/D)*(xyz/prod(xyz)^(1/D))
  fwhm <-1/resels
  
  results <-list(fwhm=fwhm,RPVImg=RPVImg,img=img)
  return(results)
}
