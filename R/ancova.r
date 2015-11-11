##############
#Prepare Data#
##############
library(ANTsR)
setwd("/Users/zach8769/Desktop/jacsobik/")

#Import MRIs into matrix
ilist<-list.files()
mask<-getMask(antsImageRead(ilist[1]))
mat<-imagesToMatrix(ilist, mask)

#Load data table
sobik<-read.table("/Users/zach8769/data/SOBIK.csv", header=TRUE,sep=",")

#Subset data to demographics and desired variable columns
vardata<-sobik[c(1:9,12)]


#Get rid of of rows missing data and create list noting which rows to keep
vardata<-na.omit(vardata)
varlist<-list(as.numeric(rownames(vardata)))

var1<-vardata[,10]
#This is my lazy way of getting the image matrix and variable 
#data table to line up. Later I'll demonstrate how to call subject
#labels from a data table so that one doesn't need to copy images
#into a single file and use more memory.
imglist<-list()
nimgs<-length(ilist)
for (i in 1:nimgs){
	imglist[i]<-smoothImage(antsImageRead(ilist[i]), 2)
	}
mat<-imageListToMatrix(imglist,mask)
for (i in varlist){
	varmat<-subset(mat[i,])
	}

##############
#ANCOVA Model#
##############

#prepare contrasts
ig<-vardata$InjuryGroup
ig<-factor(ig,levels=c("oi","mild","mod","sev"))
conmat=matrix(c(1/4,1/4,1/4,1/4,1,-1,0,0,1,0,-1,0,1,0,0,-1),ncol=4)
mymat=solve(t(conmat))
my.contrasts<-mymat[,2:4]
contrasts(ig)=my.contrasts

#Fit ancova model
voxels<-ncol(varmat)
subs<-nrow(varmat)
res<-matrix(nrow=subs, ncol=voxels)
tstat.mild<-rep(1,ncol(varmat))
tstat.mod<-rep(1,ncol(varmat))
tstat.sev<-rep(1,ncol(varmat))

progress <- txtProgressBar(min = 0, max = voxels, style = 3)
for (i in 1:voxels){
  vox<-varmat[,i]
  fit<-lm(vox~var1+ig)
  res[,i]<-residuals(fit)
  sumfit<-summary(fit)
	tstat.mild[i]<-sumfit$coefficients[3,3]
	tstat.mod[i]<-sumfit$coefficients[4,3]
	tstat.sev[i]<-sumfit$coefficients[5,3]
  setTxtProgressBar(progress, i)
  }
close(progress)

tmild<-makeImage(mask,tstat.mild)
antsImageWrite(tmild,file='/Users/zach8769/Desktop/rsobik/BASCtAggression/tmild.nii.gz')
tmod<-makeImage(mask,tstat.mod)
antsImageWrite(tmod,file='/Users/zach8769/Desktop/rsobik/BASCtAggression/tmod.nii.gz')
tsev<-makeImage(mask,tstat.sev)
antsImageWrite(tsev,file='/Users/zach8769/Desktop/rsobik/BASCtAggression/tsev.nii.gz')

rdf<-fit$df.residuals

fwhm<-estScaled.smooth(res,rdf,mask)

clusters<-rft.thresh(timg,.05,150,fwhm,mask,df,"T")

voxels<-rft.voxel(cluster$stat,fwhm,mask,df,"T")
