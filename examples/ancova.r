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
sobik<-read.table("/Users/zach8769/data/SOBIK-demographic-data.csv", header=TRUE,sep=",")

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
rdf<-fit$df.residuals

#Use the extracted statistical information for each contrast and create a statistical parametric map (SPM)
#in an antsImage format and then save it as a nifti file
tmild<-makeImage(mask,tstat.mild)
antsImageWrite(tmild,file='/Users/zach8769/Desktop/rsobik/BASCtAggression/tmild.nii.gz')
tmod<-makeImage(mask,tstat.mod)
antsImageWrite(tmod,file='/Users/zach8769/Desktop/rsobik/BASCtAggression/tmod.nii.gz')
tsev<-makeImage(mask,tstat.sev)
antsImageWrite(tsev,file='/Users/zach8769/Desktop/rsobik/BASCtAggression/tsev.nii.gz')


#Calculate the fwhm/smoothness of the image using the scaled residuals
fwhm<-estScaled.smooth(res,rdf,mask)

#####
#RFT#
#####

#So far I've created three images full of statistical information but none of it is ready to be reported.
#Because I've literally performed over a million tests a .05 probability level will give me over 50,000 
#significant voxels just by chance. I need to use a method that accounts for the many tests I perform and
#the inherent neuroanatomical relationship between proximal voxels. This is why I'll be using Random Field
#Theory (RFT) to control for multiple tests.

#Calculate an expected threshold given a certain p-value, suprathreshold cluster volume, and overall SPM size
stat<-rft.thresh(tsev,.05,150,fwhm,mask,df,"T")
#Using the previously obtained threshold find the cluster-level and voxel-level statistics
results<-rft.results(tsev,stat,150,fwhm,mask,rdf,"T")

#Extract visual information for the first cluster along the x, y, and z axis.
#For this particular fitted model only positive "T" values survived thresholding, so I'll be 
#looking at the PositiveStatistics table within results object
invisible(plot(t1,list(clust),slices=c(results$PositiveStatistics[1,5]),axis=1,overlay.color=c(red)))
invisible(plot(t1,list(clust),slices=c(results$PositiveStatistics[1,6]),axis=2,overlay.color=c(red)))
invisible(plot(t1,list(clust),slices=c(results$PositiveStatistics[1,7]),axis=3,overlay.color=c(red)))

#I like the coronal view provided in the third option so I create a jpeg and...
jpeg('/Users/zach8769/Desktop/rsobik/BASCtAggression/PCluster1.jpeg')
invisible(plot(t1,list(clust),slices=c(results$PositiveStatistics[1,7]),axis=3,overlay.color=c(red)))
dev.off()
# ...a csv of my data.
write.csv(results$PositiveStatistics, file='/Users/zach8769/Desktop/rsobik/BASCtAggresion/PStatsSev.csv')
