library(MASS)

## Simple Linear Regression

# Set up data
nsub <-nrow(imat)
nvox <-ncol(imat)

# Create Design Matrix
dm <-model.matrix(~var1-1)
dm <-cbind(dm,1)
degf <-nsub - ncol(dm)

# Linear regression through matrices
UU <-ginv(t(dm) %*% dm)
UY <-t(dm) %*% imat
B <- UU %*% UY
residuals <-imat - (dm %*% B)
rss <-colSums(res^2)
mrss <-rss
se <-sqrt(mrss *(contrast%*% UU %*% contrast))
tfield <-(contrast %*% B)/se

# Estimating the smoothness
Mmat <-colMeans(residuals)
Zmat <-matrix(nrow=nsub, ncol=nvox)
cat("Estimating fwhm/smoothing",sep="")
progress <- txtProgressBar(min = 0, max = nsub, style = 3)
for (i in 1:nsub){
	Zmat[i,]<-(res[i,]-Mmat[1])/psd
	img<-makeImage(mask,Zmat[i,])
	smooth<-est.Smooth(img,mask)
	fwhm<-fwhm+smooth[[2]]
	setTxtProgressBar(progress, i)
	}
close(progress)
fwhm2<-sqrt(4*log(2)/(fwhm/degf)

# multiple regression
dm <-model.matrix(~var1+var2)
dm <-cbind(dm,1)

# anova-one way
dm <-model.matrix(~ig-1)
dm <-cbind(dm,1)

# anova-two-way
dm <-model.matrix(~ig:im-1)
dm <-cbind(dm,1)

# ancova
dm <-model.matrix(~var1:ig-1)
dm <-cbind(dm,1)

# mancova
dm <-model.matrix(~var1:ig:im)
dm <-cbind(dm,1)

# adding controls (var2 is the control variable)
dm <-model.matrix(~ig+var2-1)
dm <-cbind(dm,1)
