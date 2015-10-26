est.smooth<-function(ilist,field){

dimx <- dim(ilist[1])[1]
dimy <- dim(ilist[1])[2]
dimz <- dim(ilist[1])[3]
n<-length(ilist)
N<-length(field)
Zxp<-0
Zyp<-0
Zzp<-0
Zxy<-0
Zyz<-0
Zxz<-0
for (i in 1:n){
	ilist[i]<-antsImageRead(ilist[i])
	}
for (i in 1:n){
	for (x in 1:(dimx)){
		for (y in 1:(dimy)){
			for (z in 1:(dimz)){
				vox<-getPixels(ilist[i],x,y,z)[1]
				Zxx<-Zxx+(((getPixels(ilist[i],x+1,y,z)[1]-vox)^2)/(N*(n-1)))
				Zyy<-Zyy+(((getPixels(ilist[i],x,y+1,z)[1]-vox)^2)/(N*(n-1)))
				Zzz<-Zzz+(((getPixels(ilist[i],x,y,z+1)[1]-vox)^2)/(N*(n-1)))
				Zxy<-Zxy+(((vox+getPixels(ilist[i],x,y+1,z)[1])*(vox+getPixels(ilist[i],x+1,y,z)[1]))/(4*N*(n-1)))
				Zxz<-Zxz+(((vox+getPixels(ilist[i],x,y,z+1)[1])*(vox+getPixels(ilist[i],x+1,y,z)[1]))/(4*N*(n-1)))
				Zyz<-Zyz+(((vox+getPixels(ilist[i],x,y+1,z)[1])*(vox+getPixels(ilist[i],x,y,z+1)[1]))/(4*N*(n-1)))		}
			}
		}
	}
fwhmx<-(4*log(2/Zxx))^(1/2)
fwhmy<-(4*log(2/Zyy))^(1/2)
fwhmz<-(4*log(2/Zzz))^(1/2)
fwhm<-c(fwhmx,fwhmy,fwhmz)
return(fwhm)
}
