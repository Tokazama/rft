est.smooth<-function(ilist,field){
n<-length(ilist)
imglist<-list()
for (i in 1:n){
	imglist[i]<-antsImageRead(ilist[i])
	}
dimx <- dim(imglist[[1]])[1]
dimy <- dim(imglist[[1]])[2]
dimz <- dim(imglist[[1]])[3]
N<-length(field)
Zxx<-0
Zyy<-0
Zzz<-0
Zxy<-0
Zyz<-0
Zxz<-0
for (i in 1:n){
	for (x in 1:(dimx)){
		for (y in 1:(dimy)){
			for (z in 1:(dimz)){
				vox<-getPixels(imglist[[i]],x,y,z)[1]
				Zxx<-Zxx+(((getPixels(imglist[[i]],x+1,y,z)[1]-vox)^2)/(N*(n-1)))
				Zyy<-Zyy+(((getPixels(imglist[[i]],x,y+1,z)[1]-vox)^2)/(N*(n-1)))
				Zzz<-Zzz+(((getPixels(imglist[[i]],x,y,z+1)[1]-vox)^2)/(N*(n-1)))
				Zxy<-Zxy+(((vox+getPixels(imglist[[i]],x,y+1,z)[1])*(vox+getPixels(imglist[[i]],x+1,y,z)[1]))/(4*N*(n-1)))
				Zxz<-Zxz+(((vox+getPixels(imglist[[i]],x,y,z+1)[1])*(vox+getPixels(imglist[[i]],x+1,y,z)[1]))/(4*N*(n-1)))
				Zyz<-Zyz+(((vox+getPixels(imglist[[i]],x,y+1,z)[1])*(vox+getPixels(imglist[[i]],x,y,z+1)[1]))/(4*N*(n-1)))		}
			}
		}
	}
fwhmx<-(4*log(2/Zxx))^(1/2)
fwhmy<-(4*log(2/Zyy))^(1/2)
fwhmz<-(4*log(2/Zzz))^(1/2)
fwhm<-c(fwhmx,fwhmy,fwhmz)
return(fwhm)
}
