est.smooth<-function(ilist,field){

dimx <- dim(antsImageRead(ilist[1]))[1]
dimy <- dim(antsImageRead(ilist[1]))[2]
dimz <- dim(antsImageRead(ilist[1]))[3]
n<-length(ilist)
N<-length(stat)
Zxp<-0
Zyp<-0
Zzp<-0
for (i in ilist){
	for (x in 1:(dimx)){
		for (y in 1:(dimy)){
			for (z in 1:(dimz)){
				vox<-getPixels(antsImageRead(ilist[i]),x,y,z)[1]
				Zxx<-Zxx+(((getPixels(antsImageRead(ilist[i]),x+1,y,z)[1]-vox)^2)/(N*(n-1)))
				Zyy<-Zyy+(((getPixels(antsImageRead(ilist[i]),x,y+1,z)[1]-vox)^2)/(N*(n-1)))
				Zzz<-Zzz+(((getPixels(antsImageRead(ilist[i]),x,y,z+1)[1]-vox)^2)/(N*(n-1)))
				Zxy<-Zxy+(((vox+getPixels(antsImageRead(ilist[i]),x,y+1,z)[1])*(vox+getPixels(antsImageRead(ilist[i]),x+1,y,z)[1]))/(4*N*(n-1)))
				Zxz<-Zxz+(((vox+getPixels(antsImageRead(ilist[i]),x,y,z+1)[1])*(vox+getPixels(antsImageRead(ilist[i]),x+1,y,z)[1]))/(4*N*(n-1)))
				Zyz<-Zyz+(((vox+getPixels(antsImageRead(ilist[i]),x,y+1,z)[1])*(vox+getPixels(antsImageRead(ilist[i]),x,y,z+1)[1]))/(4*N*(n-1)))
				}
			}
		}
	}
fwhmx<-(4*log(2/Zxx))^(1/2)
fwhmy<-(4*log(2/Zyy))^(1/2)
fwhmz<-(4*log(2/Zzz))^(1/2)
}
