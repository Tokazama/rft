estScaledSmooth <-function(res,df,mask){
  D <-mask@dimension
  n <-nrow(res)
  scale <-(1/n)*((df[1]+df[2]+1)/df[2])
  sr <-res/sqrt((res*res)/df[2])
  cat("Estimating fwhm/smoothing")
  maskar <-as.array(mask)
  voxels <-sum(maskar)
  dimx <-dim(mask)[1]
  dimy <-dim(mask)[2]
  dimz <-dim(mask)[3]
  dx2<-dy2<-dz2<-dxy<-dxz<-dyz<-array(0L,dim=c(dim(mask)+2))
  d1 <-array(0, dim=c(dim(mask)+2))
  m1 <-array(0, dim=c(dim(mask)+2))
  d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
  m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
  
  
  progress <- txtProgressBar(min = 0, max = nrow(res), style = 3)
  for (i in 1:n){
    img <-makeImage(mask,sr[,i])
    imgar <-as.array(img)

    #calculate partial derivatives x
    d2 <-array(0, dim=c(dim(mask)+2))
    m2 <-array(0, dim=c(dim(mask)+2))

    d2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
    m2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
    m3 <-(m1 + m2) == 2
    dx <-(d1 - d2)*m3

    #calculate partial derivatives y
    d2 <-array(0, dim = dim(mask) + 2)
    m2 <-array(0, dim = dim(mask) + 2)
    
    d2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
    m2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
    m3 <-(m1 + m2) == 2
    dy <-(d1 - d2)*m3

    #calculate partial derivatives z
    d2 <-array(0, dim = dim(mask) + 2)
    m2 <-array(0, dim = dim(mask) + 2)
    
    d2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
    m2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
    m3 <-(m1 + m2) == 2
    dz <-(d1 - d2)*m3
    
    dx2 <-dx2 + sum(dx*dx)/voxels
    dy2 <-dy2 + sum(dy*dy)/voxels
    dz2 <-dz2 + sum(dz*dz)/voxels
#    dxy <-dxy + sum(dx*dy)/voxels
#    dxz <-dxz + sum(dx*dz)/voxels
#    dyz <-dyz + sum(dy*dz)/voxels
    
    setTxtProgressBar(progress, i)
    }
  close(progress)
      
  # this should sum all of them into one image create variance covariance
  Vxx <-dx2 *scale
  Vyy <-dy2 *scale
  Vzz <-dz2 *scale
#  Vxy <-dxy *scale
#  Vxz <-dxz *scale
#  Vyz <-dyz *scale

  lambda <-cbind(Vxx,Vyy,Vzz)
  
  lambda <-sqrt(4*log(2)/lambda)
  fwhm <-colSums(lambda)/voxels
  
#  fwhm_img <-Vxx*Vyy*Vzz+Vxy*Vyz*Vxz*2-Vxx*Vyz*Vyz-Vxy*Vxy*Vzz-Vxz*Vyy*Vxz
  
#  fwhm_img[fwhm_img < 0] <-0
  
#  lambda <-sqrt(lambda/(4*log(2)))
#  fwhm_img <-sqrt(fwhm_img/(4*log(2))^D)
  
#  lambda <-colSums(lambda)/voxels
#  fwhm_img <-sum(fwhm_img)/voxels

#  RESEL <-fwhm_img^(1/D)*(lambda/prod(lambda)^(1/D))
#  fwhm <-1/RESEL
  fwhm
  }
