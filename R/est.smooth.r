#psdf-pooled standard deviation
#
est.smooth<-function(mat,mask,psdf){
	voxels<-ncol(mat)
	subs<-length(mat)
	Mmat<-colMeans(mat)
	dimx<-dim(mask)[1]
	dimy<-dim(mask)[2]
	dimz<-dim(mask)[3]
	for (i in 1:subs){
		Zmat[i,]<-(mat[i,]-Mmat)/psdf
		}
	lambda<-matrix(nrow=3,ncol=3)
	for (i in 1:subs){
		img<-makeImage(mask,Zmat[i,]
		for (x in 1:(dimx)){
			for (y in 1:(dimy)){
				for (z in 1:(dimz)){
					suppressWarnings(warning("index not inside the image : "[ x, y, z, ]))
					vox<-getPixels(img,x,y,z)[1]
					lambda[1,1]<-lambda[1,1]+(((getPixels(img,x+1,y,z)[1]-vox)^2)/(voxels*(subs-1)))
					lambda[2,2]<-lambda[2,2]+(((getPixels(img,x,y+1,z)[1]-vox)^2)/(voxels*(subs-1)))
					lambda[3,3]<-lambda[3,3]+(((getPixels(img,x,y,z+1)[1]-vox)^2)/(voxels*(subs-1)))
					lambda[1,2]<-lambda[1,2]+(((vox+getPixels(img,x,y+1,z)[1])*(vox+getPixels(img,x+1,y,z)[1]))/(4*voxels*(subs-1)))
					lambda[1,3]<-lambda[1,3]+(((vox+getPixels(img,x,y,z+1)[1])*(vox+getPixels(img,x+1,y,z)[1]))/(4*voxels*(subs-1)))
					lambda[3,2]<-lambda[3,2]+(((vox+getPixels(img,x,y+1,z)[1])*(vox+getPixels(img,x,y,z+1)[1]))/(4*voxels*(subs-1)))		
					}
				}
			}
		}
	lambda[2,1]<-lambda[1,2]
	lambda[3,1]<-lambda[1,3]
	lambda[2,3]<-lambda[3,2]

	fwhmx<-(4*log(2/lambda[1,1]))^(1/2)
	fwhmy<-(4*log(2/lambda[2,2]))^(1/2)
	fwhmz<-(4*log(2/lambda[3,3]))^(1/2)
	fwhm<-c(fwhmx,fwhmy,fwhmz)
##Tpdf<-probability density function of a t-distribution with v degress of freedom
	vintegral<-function(t){((((t^2)+subs-1)^2)/((v-1)*(v-2)))*((dt(t,df)^3)/(p(t)^2))}
	lambdav<-integrate(vintegral,Inf,Inf)
	lambdaz<-lambdav*lambdae
	return(fwhm,lambda)
	}
