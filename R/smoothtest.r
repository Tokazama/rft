
maskar <-as.array(mask)
voxels <-sum(maskar)
imgar <-as.array(img)
fwhm <-matrix(0L,nrow=1,ncol=3)
lambda <-matrix(0L,nrow=3,ncol=3)

dimx <-dim(mask)[1]
dimy <-dim(mask)[2]
dimz <-dim(mask)[3]
d1 <-array(0, dim = dim(mask) + 2)
m1 <-array(0, dim = dim(mask) + 2)
m2 <-array(0, dim = dim(mask) + 2)
m1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-maskar
xm2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-maskar
ym2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-maskar
zm2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-maskar
xm3 <-(m1 + xm2) == 2
ym3 <-(m1 + ym2) == 2
zm3 <-(m1 + zm2) == 2

for (i in 1:n){
#calculate partial derivatives x
d2 <-array(0, dim = dim(mask) + 2)
d1[2:(dimx+1), 2:(dimy+1), 2:(dimz+1)] <-imgar
d2[1:(dimx), 2:(dimy+1), 2:(dimz+1)] <-imgar
dx <-((d1 - d2)*xm3)
#calculate partial derivatives y
d2 <-array(0, dim = dim(img) + 2)
d2[2:(dimx+1), 1:(dimy), 2:(dimz+1)] <-imgar
dy <-((d1 - d2)*ym3)
#calculate partial derivatives z
d2 <-array(0, dim = dim(img) + 2)
d2[2:(dimx+1), 2:(dimy+1), 1:(dimz)] <-imgar
dz <-((d1 - d2)*zm3)

dxy <-(dx[1:dimx,1:dimy,1:dimz]+dx[1:dimx,2:(dimy+1),1:dimz])*ym3
dyx <-(dy[1:dimx,1:dimy,1:dimz]+dy[2:(dimx+1),1:dimy,1:dimz])*xm3
dxz <-(dx[1:dimx,1:dimy,1:dimz]+dx[1:dimx,1:dimy,2:(dimz+1)])*zm3
dzx <-(dz[1:dimx,1:dimy,1:dimz]+dz[2:(dimx+1),1:dimy,1:dimz])*xm3
dyz <-(dy[1:dimx,1:dimy,1:dimz]+dy[1:dimx,1:dimy,2:(dimz+1)])*zm3
dzy <-(dz[1:dimx,1:dimy,1:dimz]+dz[1:dimx,2:(dimy+1),1:dimz])*ym3

Vxx <-Vxx+(dx*dx)
Vyy <-Vyy+(dy*dy)
Vzz <-Vzz+(dz*dz)
Vxy <-Vxy+(dxy*dyx)
Vxz <-Vxz+(dxz*dzx)
Vyz <-Vyz+(dyz*dzy)
}
Vxx <-Vxx/(nvox*(n-1))
Vyy <-Vyy/(nvox*(n-1))
Vzz <-Vzz/(nvox*(n-1))
Vxy <-Vxy/(4*nvox*(n-1))
Vxz <-Vxz/(4*nvox*(n-1))
Vyz <-Vyz/(4*nvox*(n-1))
