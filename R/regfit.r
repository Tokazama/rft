###RFT_Test###
library(ANTsR)
sobik<-read.table('Desktop/data/SOBIK-demographic-data.csv',header=TRUE,sep=",")
mask<-antsImageRead('Desktop/7_13_template/brainmask.nii')
setwd('Desktop/jacsobik/')
ilist<-list.files()

vardata<-sobik[c(1:10,16)]

#Get rid of of rows missing data)
vardata<-na.omit(vardata)

#Create list of the rows remaining in the data table to parse out desired images from image matrix
varlist<-rownames(vardata)
varlist<-as.numeric(varlist)
varlist<-list(varlist)

#This is my lazy way of getting the image matrix and variable 
#data table to line up and smooth. I'm pretty sure there's a better ways.

imglist<-list()
nimgs<-length(ilist)
for (i in 1:nimgs){
	imglist[i]<-smoothImage(antsImageRead(ilist[i]), 2)
	}
imat<-imageListToMatrix(imglist,mask)
for (i in varlist){
	varmat<-subset(imat[i,])
	}

#Save RAM by removing big unnecessary objects
rm(imat,imglist,nimgs,sobik)
setwd("/Users/zach8769/Desktop/rsobik/TEAChCreatureCountingTotalCorrect")

#fit general linear model
subs<-nrow(varmat)
voxels<-ncol(varmat)
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
timg<-makeImage(mask,regtstat)
antsImageWrite(timg,file="Timg.nii.gz")
	
###This might be where I could start a complete rft function
#ants.rft<-function(resmat,statmat,mask,df,StatType,statdir){

res2<-colSums(res^2)
rdf<-regfit$df.residual
S2<-res2/rdf
psd<-sqrt(S2)
Mmat<-colMeans(res)
Zmat<-res
fwhm<-matrix(0L,nrow=1,ncol=3)
subs<-nrow(res)

cat("Estimating fwhm/smoothing")
progress <- txtProgressBar(min = 0, max = subs, style = 3)
for (i in 1:subs){
  Zmat[i,]<-(res[i,]-Mmat[1])/psd
  img<-makeImage(mask,Zmat[i,])
	smooth<-est.smooth(img,mask,1,1,1)
	fwhm<-fwhm+smooth[[2]]
	setTxtProgressBar(progress, i)
	}
close(progress)
fwhm2<-sqrt(4*log(2)/(fwhm/(subs-1)))

cat("Determing clustering threshold using RFT")
thresh<-rft.thresh(pval,ka,fwhm,mask,df,fieldtype)

posclust<-image2ClusterImages(timg,150,thresh,Inf)
postable<-matrix(ncol=5)
colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
for (i in 1:length(posclust)){
	cat("Determing positive cluster level statistics:",i,sep=" ")
	cmask<-getMask(posclust[[i]])
	cvoxs<-sum(as.array(cmask))
	pclust<-rft.pcluster(posclust[[i]],mask,fwhm,thresh,df,fieldtype)
	peak<-max(posclust[[i]])
	loc<-labelImageCentroids(posclust[[i]])
	resel <-ants.resel(mask,fwhm)
	ec<-ants.ec(stat,fieldtype,df)
	pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
	postable[i,]<-c(cvox, pclust, peak, pval, loc)
	clustername<-paste("P-Cluster:",i,sep="")
	rownames(postable[i,])<-c(clustername)
	}

negclust<-image2ClusterImages(timg,150,-Inf,-thresh)
negtable<-matrix(ncol=5)
colnames(postable)<-c("Voxels", "Cluster-Probability", "Peak-Height", "Voxel-Probability", "Coordinates")
for (i in 1:length(posclust)){
	cat("Determing negative cluster level statistics:",i,sep=" ")
	cmask<-getMask(negclust[[i]])
	cvoxs<-sum(as.array(cmask))
	pclust<-rft.pcluster(posclust[[i]],mask,fwhm,thresh,df,fieldtype)
	peak<-max(posclust[[i]])
	loc<-labelImageCentroids(posclust[[i]])
	resel <-ants.resel(mask,fwhm)
	ec<-ants.ec(stat,fieldtype,df)
	pval<-(resel[1]*ec[1])+(resel[2]*ec[2])+(resel[3]*ec[3])+(resel[4]*ec[4])
	negtable[i,]<-c(cvox, pclust, peak, pval, loc)
	clustername<-paste("N-Cluster:",i,sep="")
	rownames(postable[i,])<-c(clustername)
	}
cluster.table<-rbind(postable,negtable)

	
}
	



###Try to imitate cluster level statistic

###Try to imitate voxel level statistic

###Put together all 


 
  
  
