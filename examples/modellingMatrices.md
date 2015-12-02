Random Field Theory
===================
Once jacobian images are obtained from a morphometry pipeline there are a variety of statistical models that can be used to evaluate these images. Some pipelines warp a template atlas along with the template brain into subject space. Volumes can be obtained from each region of interest (ROI) in the warped atlases. Another common approach is to perform a statistical test at each individual voxel. Modern magnetic resonance images (MRIs) typically contain over a million voxels in a single scan. From a statistical standpoint, this poses a challenge. The probability of a false-positive necesitates an alpha far below the traditional .05. If one were to perform one million statistical tests with an alpha of .05 and using traditional bonferroni correction the adjusted alpha would be .05e-6. 

There are a variety of ways 
Random field theory (RFT) allows us to take statistical results run on voxel-based morphometry (VBM) data and more realistically threshold 

Simple Linear Regression
------------------------

### Set up data

```
library(ANTsR)
library(MASS)
nsub <-nrow(imat)
nvox <-ncol(imat)
```
### Create Design Matrix
```
dm <-model.matrix(~var1-1)
dm <-cbind(dm,1)
degf <-nsub - ncol(dm)
```

### Calculating T-Scores
```
UU <-ginv(t(dm) %*% dm)
UY <-t(dm) %*% imat
B <- UU %*% UY
residuals <-imat - (dm %*% B)
rss <-colSums(res^2)
mrss <-rss
se <-sqrt(mrss *(contrast%*% UU %*% contrast))
tfield <-(contrast %*% B)/se
```

### Estimating the smoothness
```
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
```

Multiple Regression
-------------------
```
dm <-model.matrix(~var1+var2)
dm <-cbind(dm,1)
```

Anova (one-way)
---------------

```
dm <-model.matrix(~ig-1)
dm <-cbind(dm,1)
```

Anova (two-way)
---------------
```
dm <-model.matrix(~ig:im-1)
dm <-cbind(dm,1)
```

Ancova
------
```
dm <-model.matrix(~var1:ig-1)
dm <-cbind(dm,1)
```

Mancova
-------
```
dm <-model.matrix(~var1:ig:im)
dm <-cbind(dm,1)
```

Adding Controls
---------------
(var2 is the control variable)
```
dm <-model.matrix(~ig+var2-1)
dm <-cbind(dm,1)
```
