library(ANTsR)
setwd("/fslhome/zach8769/compute/jacsobik/")
ilist<-list.files()
mask<-getMask(antsImageRead(ilist[1]))
mat<-imagesToMatrix(ilist, mask)

#Load data table
sobik<-read.table("Desktop/data/sobik.csv", header=TRUE,sep=",")
#Subset data to demographics and desired variable columns
vardata<-sobik[c(1:9,12)]
#Get rid of of rows missing data)
vardata<-na.omit(vardata)
#Create list of the rows remaining in the data table
varlist<-rownames(vardata)
varlist<-as.numeric(varlist)
varlist<-list(varlist)
#
for(i in varlist){
  varmat<-subset(mat[i,])
  }
#prepare contrasts
ig<-vardata$InjuryGroup
ig<-factor(ig,levels=c("oi","mild","mod","sev"))
conmat=matrix(c(1/4,1/4,1/4,1/4,1,-1,0,0,1,0,-1,0,1,0,0,-1),ncol=4)
mymat=solve(t(conmat))
my.contrasts<-mymat[,2:4]
contrasts(ig)=my.contrasts
var1<-vardata[,10]
fit<-lm(varmat~var1+ig)
save(fit,file="fit")
