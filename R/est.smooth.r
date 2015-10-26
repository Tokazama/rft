est.smooth<-function(lmfit,mat,mask){
	res<-residuals(lmfit)
	res2<-colSums(res^2)
	rdf<-lmfit$df.residual
	S2<-res2/redf
	S<-sqrt(S2)
	voxels<-ncol(mat)
	subs<-length(mat)
	Mmat<-matrix(nrow=1,ncol=voxels)
	for (i in 1:voxels){
		Mmat[1,i]<-colMeans(mat[,i])
		}
	Zmat<-(mat-Mmat)/S
	dimx<-dim(mask)[1]
	dimy<-dim(mask)[2]
	dimz<-dim(mask)[3]
	for (i in 1:subs){
		Zmat[,i]<-(mat[,i]-Mmat[,1])/(
		}
	Zxx<-0
	Zyy<-0
	Zzz<-0
	Zxy<-0
	Zyz<-0
	Zxz<-0
	for (i in 1:subs){
		img<-makeImage(mask,Zmat[i,]
		for (x in 1:(dimx)){
			for (y in 1:(dimy)){
				for (z in 1:(dimz)){
					vox<-getPixels(img,x,y,z)[1]
					Zxx<-Zxx+(((getPixels(img,x+1,y,z)[1]-vox)^2)/(voxels*(subs-1)))
					Zyy<-Zyy+(((getPixels(img,x,y+1,z)[1]-vox)^2)/(voxels*(subs-1)))
					Zzz<-Zzz+(((getPixels(img,x,y,z+1)[1]-vox)^2)/(voxels*(subs-1)))
					Zxy<-Zxy+(((vox+getPixels(img,x,y+1,z)[1])*(vox+getPixels(img,x+1,y,z)[1]))/(4*voxels*(subs-1)))
					Zxz<-Zxz+(((vox+getPixels(img,x,y,z+1)[1])*(vox+getPixels(img,x+1,y,z)[1]))/(4*voxels*(subs-1)))
					Zyz<-Zyz+(((vox+getPixels(img,x,y+1,z)[1])*(vox+getPixels(img,x,y,z+1)[1]))/(4*voxels*(subs-1)))		
					}
				}
			}
		}
	lambda<-dim(nrow=3,ncol=3)
	lambda[1,1]<-(4*log(2/Zxx))^(1/2)
	lambda[2,2]<-(4*log(2/Zyy))^(1/2)
	lambda[3,3]<-(4*log(2/Zzz))^(1/2)
	lambda[1,2]<-(4*log(2/Zxy))^(1/2)
	lambda[2,1]<-(4*log(2/Zxy))^(1/2)
	lambda[1,3]<-(4*log(2/Zxz))^(1/2)
	lambda[3,1]<-(4*log(2/Zxz))^(1/2)
	lambda[2,3]<-(4*log(2/Zyz))^(1/2)
	lambda[3,2]<-(4*log(2/Zyz))^(1/2)
	fwhm<-c(fwhmx,fwhmy,fwhmz)
	lambda<-dim(nrow=3,ncol=3)
##Tpdf<-probability density function of a t-distribution with v degress of freedom
	vintegral<-function(t){((((t^2)+subs-1)^2)/((v-1)*(v-2)))*((Tpdf(t)^3)/(p(t)^2))}
	lambdav<-integrate(vintegral,Inf,Inf)
	lambdaz<-lambdav*lambdae
	return(fwhm)
	}
