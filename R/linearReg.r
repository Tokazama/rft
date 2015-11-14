###RFT_Test###
library(ANTsR)
sobik<-read.table('Desktop/data/SOBIK-demographic-data.csv',header=TRUE,sep=",")
mask<-antsImageRead('Desktop/7_13_template/brainmask.nii')

vardata<-sobik[c(1:10,16)]

#Get rid of of rows missing data)
vardata<-na.omit(vardata)

#The first column of the sobik data table contains the respective folder names for each subject.
#This allows me to put together a fairly simple process in which I can read in any image type from
#the subject folders into the same order as my variables without having to worry about messing
#with the order. I also smooth the images before I convert them into a matrix.
subs<-sobik[,1]
paths<-paste("/Users/zach8769/Desktop/sobik/",subs,"/warp/mod.nii.gz",sep="")
ilist<-imageFileNames2ImageList(paths)
imglist<-list()
for (i in 1:length(ilist)){
	imglist[i]<-smoothImage(ilist[[i]],2)
	}
imat<-imageListToMatrix(ilist,mask)

#fit general linear model
subs<-nrow(imat)
voxels<-ncol(imat)
regtstat<-matrix(nrow=1,ncol=voxels)
res<-matrix(0L,nrow=subs,ncol=voxels)
var1<-vardata[,10]

progress <- txtProgressBar(min = 0, max = voxels, style = 3)
for (i in 1:voxels){
	vox<-varmat[,i]
	regfit<-lm(vox~var1)
	##Extract statistical values
	res[,i]<-residuals(regfit)
	regsum<-summary(regfit)
	regtstat[,i]<-regsum$coefficients[2,3]
	setTxtProgressBar(progress, i)
	}
	close(progress)
rdf<-regfit$df.residuals
timg<-makeImage(mask,regtstat)
antsImageWrite(timg,file="treg.nii.gz")

fwhm<-estScaled.smooth(res,rdf,mask)

stat<-rft.thresh(tsev,.05,150,fwhm,mask,df,"T")
results<-rft.results(tsev,stat,150,fwhm,mask,rdf,"T")
