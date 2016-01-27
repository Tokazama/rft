'# Nonstationary cluster-size inference with random field and permutation methods
'# Hayasaka et al. (2004)
'# A Three-Dimensional Statistical Analysis for CBF Activation Studies in Human Brain
'# Worsley et al., (1992)
'# Robust Smoothness Estimation in Statistical Parametric Maps Using Standardized Residuals from the General Linear Model
'# Kiebel et al., (1999)
estSmooth <-function(x,mask,df,sample){
  if (class(x)=="antsImage"){
    classval <-1
  }else if(class(x)=="numeric"){
    classval <-2
    x <-matrix(x,nrow=1)
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
    scale <-(1/n)*(nfull/(df[2]))
    # standardized residuals
    ux <-x/sqrt(colSums(x*x)/n)
  }
  maskar <-as.array(mask)
  nvox <-sum(maskar)
  fwhm <-matrix(0L,nrow=1,ncol=3)
  lambda <-matrix(0L,nrow=3,ncol=3)
  
  # used in find partial derivatives (taken out of loop to save time)
  dimx <-dim(mask)[1]
  dimy <-dim(mask)[2]
  dimz <-dim(mask)[3]
  d1<-dx2<-dy2<-dz2<-m1<-xm2<-ym2<-zm2<-Vxx<-Vyy<-Vzz<-Vxy<-Vxz<-Vyz<-array(0, dim = dim(mask) + 2)
  maskar <-as.array(mask)
  m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
  xm2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
  ym2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
  zm2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
  xm3 <-(m1 + xm2) == 2
  ym3 <-(m1 + ym2) == 2
  zm3 <-(m1 + zm2) == 2
  rm(xm2,ym2,zm2)
  progress <-txtProgressBar(min=0, max=n, style=3)
  for (i in 1:n){
    if (classval > 1){
      imgar <-as.array(makeImage(mask,x[i,]))
    }
    #calculate partial derivatives x
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
    dx2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
    dx <-(d1 - d2)*xm3
    
    x.v0 <-d1^2
    x.v1 <-((m1-xm2)*xm3)^2
    xx <- -(1/(4*log(1-x.v1/(2*x.v0))))
    
    #calculate partial derivatives y
    dy2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
    dy <-(d1 - d2)*ym3
    
    y.v0 <-d1^2
    y.v1 <-((m1-ym2)*ym3)^2
    yy <- -(1/(4*log(1-y.v1/(2*y.v0))))
    
    #calculate partial derivatives z
    dz2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
    dz <-(d1 - d2)*zm3
    
    z.v0 <-d1^2
    z.v1 <-((m1-zm2)*zm3)^2
    zz <- -(1/(4*log(1-z.v1/(2*z.v0))))
    
    #calculate partial derivatives yx
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-dy[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)]
    dx2[1:dimx,2:(dimy+1),2:(dimz+1)] <-dy[1:dimx,2:(dimy+1),2:(dimz+1)]
    dyx <-(d1+dx2)*xm3
    
    #calculate partial derivatives zx
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-dz[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)]
    dx2[1:dimx,2:(dimy+1),2:(dimz+1)] <-dz[1:dimx,2:(dimy+1),2:(dimz+1)]
    dzx <-(d1+dx2)*xm3
    
    #calculate partial derivatives zy
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-dz[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)]
    dy2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-dz[2:(dimx+1),1:dimy,2:(dimz+1)]
    dzy <-(d1+dy2)*ym3
    
    #calculate partial derivatives xy
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-dx[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)]
    dy2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-dx[2:(dimx+1),1:dimy,2:(dimz+1)]
    dxy <-(d1+dy2)*ym3
    
    #calculate partial derivatives yz
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-dy[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)]
    dz2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-dy[2:(dimx+1),2:(dimy+1),1:dimz]
    dyz <-(d1+dz2)*zm3
    
    #calculate partial derivatives xz
    d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-dx[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)]
    dz2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-dx[2:(dimx+1),2:(dimy+1),1:dimz]
    dxz <-(d1+dz2)*zm3
    
    Vxx <-Vxx+(dx*dx)
    Vyy <-Vyy+(dy*dy)
    Vzz <-Vzz+(dz*dz)
    Vxy <-Vxy+(dxy*dyx)
    Vxz <-Vxz+(dxz*dzx)
    Vyz <-Vyz+(dyz*dzy)
    setTxtProgressBar(progress,i)
  }
  close(progress)
  rm(dxy,dxz,dyx,dzx,dyz,dzy,dx,dy,dz,xm3,ym3,zm3)
  Vxx <-Vxx/(n-1)
  Vyy <-Vyy/(n-1)
  Vzz <-Vzz/(n-1)
  Vxy <-Vxy/(4*(n-1))
  Vxz <-Vxz/(4*(n-1))
  Vyz <-Vyz/(4*(n-1))

  
  resel.img <-(Vxx*Vyy*Vzz)+(Vxy*Vyz*Vxz*2)-(Vxx*Vyz*Vyz)-(Vxy*Vxy*Vzz)-(Vxz*Vyy*Vxz)
  resel.img[resel.img < 0] <-0
  resel.img <-sqrt(resel.img/(4*log(2))^D)
  resel.img <-sum(resel.img)/nvox
  
  fwhm <-cbind(Vxx,Vyy,Vzz)
  fwhm <-sqrt(fwhm/(4*log(2)))
  fwhm <-colSums(fwhm)/nvox
  
  resels <-resel.img^(1/D)*(fwhm/prod(fwhm)^(1/D))
  fwhm <-1/resels
  names(fwhm) <-c("x","y","z")
  return(fwhm)
  } 
