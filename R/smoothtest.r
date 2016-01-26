

dimx <-dim(img)[1]
dimy <-dim(img)[2]
dimz <-dim(img)[3]
d1 <-array(0, dim = dim(img) + 2)
d2 <-array(0, dim = dim(img) + 2)
m1 <-array(0, dim = dim(img) + 2)
m2 <-array(0, dim = dim(img) + 2)
maskar <-as.array(mask)
voxels <-sum(maskar)
imgar <-as.array(img)
fwhm <-matrix(0L,nrow=1,ncol=3)
lambda <-matrix(0L,nrow=3,ncol=3)

#calculate partial derivatives x
d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
d2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
m2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
m3 <- (m1 + m2) == 2
dx <-((d1 - d2)*m3)


#calculate partial derivatives y
d2 <-array(0, dim = dim(img) + 2)
m2 <-array(0, dim = dim(img) + 2)

d2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
m2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
m3 <-(m1 + m2) == 2
dy <-((d1 - d2)*m3)


#calculate partial derivatives z
d2 <-array(0, dim = dim(img) + 2)
m2 <-array(0, dim = dim(img) + 2)

d2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
m2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
m3 <-(m1 + m2) == 2
dz <-((d1 - d2)*m3)

Vxx <-(dx*dx)/(nvox*df[2])
Vyy <-(dy*dy)/(nvox*df[2])
Vzz <-(dz*dz)/(nvox*df[2])
Vxy <-((dx[1:dimx,1:dimy,1:dimz]+dx[1:dimx,2:(dimy+1),1:dimz])*(dy[1:dimx,1:dimy,1:dimz]+dy[2:(dimx+1),1:dimy,1:dimz]))/(4*nvox*df[2])
Vxz <-((dx[1:dimx,1:dimy,1:dimz]+dx[1:dimx,1:dimy,2:(dimz+1)])*(dz[1:dimx,1:dimy,1:dimz]+dz[2:(dimx+1),1:dimy,1:dimz]))/(4*nvox*df[2])
Vyz <-((dy[1:dimx,1:dimy,1:dimz]+dy[1:dimx,1:dimy,2:(dimz+1)])*(dz[1:dimx,1:dimy,1:dimz]+dz[1:dimx,2:(dimy+1),1:dimz]))/(4*nvox*df[2])

