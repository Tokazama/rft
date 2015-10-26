library(ANTsR)
setwd("/fslhome/zach8769/compute/rsobik/WASIMatrixReasoning/")
load("/fslhome/zach8769/compute/rsobik/WASIMatrixReasoning/vardata")
load("/fslhome/zach8769/compute/rsobik/WASIMatrixReasoning/varmat")
mask<-getMask(antsImageRead("/fslhome/zach8769/compute/7_13_template/brainmask.nii.gz"))
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
