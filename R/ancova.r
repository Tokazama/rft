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
sobik<-read.table("/path/to/directory/sobik.csv", header=TRUE,sep=",")

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
for(i in varlist){
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
tstat<-matrix(nrow=1, ncol=voxels)
res<-matrix(nrow=subs, ncol=voxels)

progress <- txtProgressBar(min = 0, max = voxels, style = 3)
for (i in 1:voxels){
  vox<-varmat[,i]
  fit<-lm(vox~var1+ig)
  res[,i]<-residuals(fit)
  regsum<-summary(fit)
  tstat[,i]<-regsum$coefficients[3,3]
  setTxtProgressBar(progress, i)
  }
close(progress)

rdf<-fit$df.residuals

fwhm<-estPooled.smooth(res,rdf,mask)

clusters<-rft.thresh(timg,.05,150,fwhm,mask,df,"T")

voxels<-rft.voxel(cluster$stat,"T",df,fwhm,mask)
