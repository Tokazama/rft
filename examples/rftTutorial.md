# Random Field Theory

Once jacobian images are obtained from a morphometry pipeline there are a variety of statistical models that can be used to evaluate these images. Some pipelines warp a template atlas along with the template brain into subject space. Volumes can be obtained from each region of interest (ROI) in the warped atlases. Another common approach is to perform a statistical test at each individual voxel. The resulting image is a statistical parametric map (SPM). Modern magnetic resonance images (MRIs) typically contain over a million voxels in a single scan. From a statistical standpoint, this poses a challenge. The probability of a false-positive necesitates an alpha far below the traditional .05. If one were to perform one million statistical tests with an alpha of .05 and using traditional bonferroni correction the adjusted alpha would be .05e-6. 

The greatest issue with this approach is that it assumes each voxel to be completely independent from neighboring voxels. In order to correct for this error one may use a priori masking to limit the region of tested voxels. This naturally decreases the undesirable effect of most statistical corrections, such as bonferroni. However, such techniques still rely upon the false assumption that each voxel is completely independent and preclude the possibility of investigative whole-brain analyses. 

Random field theory (RFT) allows us to take SPM images and consider each voxel as a component of its region rather than completely independent. Although I will focus on evaluating SPMs obtained from VBM methods herein, RFT principles can be applied to a variety of other image modalities (i.e. fMRI, PET, etc). The difficulty in utilizing RFT is that each component of the analysis on the warped images is important to the final results.

1. Obtaining residuals from the statistical model
2. Using image residuals in order to estimate smoothness
3. Calculating image resels
4. Clustering the SPM image by a chosen threshold
5. Determining cluster-level statistics
6. Determining voxel-level statistics

Because this process can become overwhelming to anyone approaching it for the first time with a functional adrenal gland, I will go through each step as thoroughly as I can. Feel free to skip ahead as you come across fairly rudimentary material.

## Fitting your data to a statistical model

The following will highlight the initial steps in obtaining residuals and a SPM images using statistical models. The important thing to realize is that none of these methods are absolute. They should merely be used as an introduction to how you can approach this step. The benefit to performing a RFT based VBM analysis in R is that you get to decide what assumptions to make. For example, additional steps of scaling images or incorporating pooled deviations between conditions and image voxels may be beneficial.

### Using `lm()`

For seasoned R users many already existing functions are sufficient to extract pertinent information. 

Here we'll set up our data for this tutorial
```
subs <-nrow(imat)
nvox <-ncol(imat)

var1 <-vardata[,10]
var2 <-
cond1 <-vardata
cond2 <-vardata
```

The following performs a statistical test utilizing the `lm()` function to acquire residuals and a SPM of T-scores. It simply demonstrates how R can be utilized as is to obtain the information we need. Before pursuing your own analysis it's encouraged that you explore other methods of obtaining the SPM and residuals that better fits the goal of your research.

```
spm <-matrix(nrow=1,ncol=voxels)
res <-matrix(0L,nrow=subs,ncol=voxels)
progress <- txtProgressBar(min = 0, max = voxels, style = 3)
for (i in 1:nvox){
	vox<-imat[,i]
	regfit<-lm(vox~var1)
	res[,i]<-residuals(regfit)
	regsum<-summary(regfit)
	spm[,i]<-regsum$coefficients[2,3]
	setTxtProgressBar(progress, i)
	}
	close(progress)
timg<-makeImage(mask, spm)
```

### Matrix Based Regression

Traditionally, software that has utilized RFT uses design matrices to quickly organize and compute results. In this subsection I will demonstrate how to perform a simple linear regression (similar to what was just performed) and how to create design matrices for other statistical models.

Here we see how we can create a disgn matrix and extract the degrees of freedom we will need in later steps.

```
dm <-model.matrix(~var1-1)
dm <-cbind(dm,1)
degf <-nsub - ncol(dm)
```

Once a design matrix is available we need to make a contrast matrix. This is similar to other methods that require a contrast matrix in R. This is fairly simple for a linear regression.

```
conmat <-c(1,0)
```

The following demonstrates the matrix manipulation necessary to extract the SPM and residuals.

```
UU <-ginv(t(dm) %*% dm)
UY <-t(dm) %*% imat
B <- UU %*% UY
res <-imat - (dm %*% B)
rss <-colSums(res^2)
mrss <-rss
se <-sqrt(mrss *(conmat%*% UU %*% conmat))
spm <-(conmat %*% B)/se
```
We can similrly create design matrices for a variety of situations

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
