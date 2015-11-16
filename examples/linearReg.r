###################
#Linear Regression#
###################

library(ANTsR)
sobik<-read.table('Desktop/data/SOBIK-demographic-data.csv',header=TRUE,sep=",")
mask<-antsImageRead('Desktop/7_13_template/brainmask.nii')


#Subset the data to TBI and columns of relevent variables and get rid of rows with missing data
vardata<-sobik[sobik$InjuryGroup=="sev"|sobik$InjuryGroup=="mod"|sobik$InjuryGroup=="mild",]
vardata<-vardata[c(1:12,13)]
vardata<-na.omit(vardata)

#Alternatively, I can subset my data to just orthopedic injury (oi)
vardata<-sobik[sobik$InjuryGroup=="oi",]
vardata<-vardata[c(1:12,13)]
vardata<-na.omit(vardata)


#The first column of the sobik data table contains the respective folder names for each subject.
#This allows me to put together a fairly simple process in which I can read in any image type from
#the subject folders into the same order as my variables without having to worry about messing
#with the order. I also smooth the images before I convert them into a matrix.
mask<-antsImageRead('Desktop/7_13_template/brainmask.nii')
subs<-vardata[,1]
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
	vox<-imat[,i]
	regfit<-lm(vox~var1)
	res[,i]<-residuals(regfit)
	regsum<-summary(regfit)
	regtstat[,i]<-regsum$coefficients[2,3]
	setTxtProgressBar(progress, i)
	}
	close(progress)
rdf<-regfit$df.residual
timg<-makeImage(mask,regtstat)
antsImageWrite(timg,file="Desktop/Ttbi.nii.gz")

fwhm<-estScaled.smooth(res,rdf,mask)
stat<-rft.thresh(timg,.05,150,fwhm,mask,rdf,"T")
results<-rft.results(timg,stat,100,fwhm,rdf,"T")

antsImageWrite(results$Cluster1,file='Desktop/tertiles/tbicluster1.nii.gz')
antsImageWrite(results$Cluster2,file='Desktop/tertiles/tbicluster2.nii.gz')
write.csv(results$stats,file='Desktop/tertiles/tbi.csv')

jpeg('Desktop/tertiles/xtbiclus1.jpeg')
invisible(plot(t1,list(results$Cluster1),slices=c(results$PositiveStatistics[1,5]),axis=1,overlay.color=c(red)))
dev.off()
