library(ANTsR)
setwd("/path/to/warped/images/")

#Import MRIs into matrix
ilist<-list.files()
mask<-getMask(antsImageRead(ilist[1]))
mat<-imagesToMatrix(ilist, mask)

#Load data table
sobik<-read.table("/path/to/directory/sobik.csv", header=TRUE,sep=",")

#Subset data to demographics and desired variable columns
vardata<-sobik[c(1:9,12)]
var1<-vardata[,10]

#Get rid of of rows missing data)
vardata<-na.omit(vardata)

#Create list of the rows remaining in the data table
varlist<-rownames(vardata)
varlist<-as.numeric(varlist)
varlist<-list(varlist)

#This is my lazy way of getting the image matrix and variable 
#data table to line up. I'm pretty sure there's a better ways.
for(i in varlist){
  varmat<-subset(mat[i,])
  }

###################
#Linear Regression#
###################

regfit<-lm(varmat~var1)

##Extract statistical values
regsum<-summary(regfit)
regpval<-regsum$coefficients[3,4]
regtstat<-regsum$coefficients[3,3]

#Convert to statistical images to view
i.regpval<-makeImage(mask,regpval)
antsImageWrite(i.regpval,file="/path/to/directory/regpval.nii.gz")
i.regtstat<-makeImage(mask,regtstat)
antsImageWrite(i.regtstat,file="/path/to/directory/regtstat.nii.gz")

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
a.fit<-lm(varmat~var1+ig)
a.sum<-summary(a.fit)

#extract pvalues
pval.mild<-sumfit$coefficients[3,4]
pval.mod<-sumfit$coefficients[4,4]
pval.sev<-sumfit$coefficients[5,4]
ip.mild<-makeImage(mask,pval.mild)
antsImageWrite(ip.mild,file="/path/to/directory/pmild.nii.gz")
ip.mod<-makeImage(mask,pval.mod)
antsImageWrite(ip.mod,file="/path/to/directory/pmod.nii.gz")
ip.sev<-makeImage(mask,pval.sev)
antsImageWrite(ip.sev,file="/path/to/directory/psev.nii.gz")

#extract tstat
tstat.mild<-sumfit$coefficients[3,3]
tstat.mod<-sumfit$coefficients[4,3]
tstat.sev<-sumfit$coefficients[5,3]
it.mild<-makeImage(mask,tstat.mild)
antsImageWrite(it.mild,file="/path/to/directory/tmild.nii.gz")
it.mod<-makeImage(mask,tstat.mod)
antsImageWrite(it.mod,file="/path/to/directory/tmod.nii.gz")
it.sev<-makeImage(mask,tstat.sev)
antsImageWrite(it.sev,file="/path/to/directory/tsev.nii.gz")
