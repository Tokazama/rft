root <-"/Users/zach8769/Desktop/"
ptbp <-read.csv(paste(root,"/ptbp/data/ptbp_summary_demographics.csv",sep=""),header=TRUE,sep=",")
ptbp <-ptbp[c(1:5)]
ptbp <-na.omit(ptbp)
rows.ptbp.train <-sample(nrow(ptbp), round(.5*nrow(ptbp)), replace=F)
ptbp.train <-ptbp[rows.ptbp.train,]
age<-ptbp.train$AgeAtScan
mask <-antsImageRead(paste(root,"/ptbp/ptbp_template/PTBP_T1_BrainCerebellumMask.nii.gz",sep=""))
sub <-ptbp.train[,1]
date <-ptbp.train[,2]
paths.ptbp.train <-paste(root,"ptbp/rft/warp/", sub, "_", date, ".nii.gz",sep="")
img.list <-list()
x <-cbind(1L,age)
conmat <-matrix(c(0,1),nrow=1)
colnames(x) <-c("Intercept","Age")
colnames(conmat) <-c("Intercept","Age")

for (i in paths.ptbp.train){
	img.list[i] <-smoothImage(antsImageRead(i), 1)
	}
imat1 <-imageListToMatrix(img.list, mask)
fit1 <-rft.lm(x,imat1,conmat,mask)
fwhm1 <-fit1$fwhm
# fwhm=c(4.755874,4.439168,4.309534)
# resels=c(1 107.3829,2496.4898,14133.7936)


for (i in paths.ptbp.train){
	img.list[i] <-smoothImage(antsImageRead(i), 1.5)
	}
imat2 <-imageListToMatrix(img.list, mask)
fit2 <-rft.lm(x,imat2,conmat,mask)
fwhm2 <-fit2$fwhm

# fwhm=c(4.841590,4.543965,4.438124)
# resels=c(1,104.909.,2381.9027,13170.3890)

for (i in paths.ptbp.train){
	img.list[i] <-smoothImage(antsImageRead(i), 2)
	}
imat3 <-imageListToMatrix(img.list, mask)
fit3 <-rft.lm(x,imat3,conmat,mask)

# fwhm=c(4.977359, 4.677361, 4.586128)
# resels=c(1.0000,   101.8433,  2244.1573, 12044.1197)

## Single Image ##

smooth.table <-matrix(ncol=5)
colnames(smooth.table) <-c("Sigma","Kernel_Width","xFWHM","yFWHM","zFWHM")
sig.seq <-seq(1,10,by=.5)
kern.seq <-seq(10,140, by=10)
num <-0
for (sigma in sig.seq){
	sig.name <-paste("Sigma-",sigma,sep="")
	for (kernel in kern.seq){
		simg <-smoothImage(img,sigma,FWHM=TRUE,max_kernel_width=kernel)
		fwhm <-estSmooth(simg, mask)
		myrow <-cbind(sigma,kernel,fwhm[1],fwhm[2],fwhm[3])
		smooth.table <-rbind(smooth.table,myrow)
		num <-num+1
		cat(num/266)
		}
	}


