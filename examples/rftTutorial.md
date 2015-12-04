# Random Field Theory

Once jacobian images are obtained from a morphometry pipeline there are a variety of statistical models that can be used to evaluate these images. Some pipelines warp a template atlas along with the template brain into subject space. Volumes can be obtained from each region of interest (ROI) in the warped atlases. Another common approach is to perform a statistical test at each individual voxel. The resulting image is a statistical parametric map (SPM). Modern magnetic resonance images (MRIs) typically contain over a million voxels in a single scan. From a statistical standpoint, this poses a challenge. The probability of a false-positive necesitates an alpha far below the traditional .05. If one were to perform one million statistical tests with an alpha of .05 and using traditional bonferroni correction the adjusted alpha would be .05e-6. 

The greatest issue with this approach is that it assumes each voxel to be completely independent from neighboring voxels. In order to correct for this error one may use a priori masking to limit the region of tested voxels. This naturally decreases the undesirable effect of most statistical corrections, such as bonferroni. However, such techniques still rely upon the false assumption that each voxel is completely independent and preclude the possibility of investigative whole-brain analyses. 

Random field theory (RFT) allows us to take SPM images and consider each voxel as a component of its region rather than completely independent. Although I will focus on evaluating SPMs obtained from VBM methods herein, RFT principles can be applied to a variety of other image modalities (i.e. fMRI, PET, etc). The difficulty in utilizing RFT is that each component of the analysis on the warped images is important to the final results.

1. Obtaining residuals from the statistical model
2. Using image residuals in order to estimate smoothness
3. Clustering the SPM image by a chosen threshold
4. Calculating image resels
5. Determining cluster-level statistics
6. Determining voxel-level statistics

Because this process can become overwhelming to anyone approaching it for the first time with a functional adrenal gland, I will go through each step as thoroughly as I can. Feel free to skip ahead as you come across fairly rudimentary material.

The examples below will be ran on the data from the sobik data set in the rft repository.


## Step 1: Fitting your data to a statistical model and extracting residuals

The following will highlight the initial steps in obtaining residuals and a SPM images using statistical models. The important thing to realize is that none of these methods are absolute. They should merely be used as an introduction to how you can approach this step. The benefit to performing a RFT based VBM analysis in R is that you get to decide what assumptions to make. For example, additional steps of scaling images or incorporating pooled deviations between conditions and image voxels may be beneficial. Before performing your own analysis it's encouraged that you explore other methods of obtaining the SPM and residuals that may better fit your research goals.

### Using `lm()`

For seasoned R users many already existing functions are sufficient to extract pertinent information.


Here we'll set up our csv data in and organize it accordingly.
```
library(MASS)
library(ANTsR)
rootdir <-'/Users/zach8769/Desktop/sobik/'

setwd(rootdir)
oasis <-read.csv('oasis.csv')
oasis[,1] ## I use this to see where the first row of the data will be using is (scan 420)
oasis <-oasis[381:416,]

# We're only concerned in 4 variables and subjects without missing data, so we subset our data table.
oasis <-oasis[c(1:3,5,8,9)]
oasis <-na.omit(oasis)
# In theis example we will be treating clinical dementia ratings as factors
oasis$CDR[oasis$CDR==0] <-"notdem"  # notdemented
oasis$CDR[oasis$CDR==.5] <-"verymild" # very mild dementia
oasis$CDR[oasis$CDR==1] <-"mild"
oasis$CDR <-as.factor(oasis$CDR)
cdr <-oasis$CDR 
age <-oasis$Age
mmse <-oasis$MMSE #mini mental state exam
sex <-oasis$M.F

```

The following performs a statistical test utilizing the `lm()` function to acquire residuals and a SPM of T-scores. It simply demonstrates how R can be utilized as is to obtain the information we need. 

```
spm <-matrix(nrow=1,ncol=voxels)
res <-matrix(0L,nrow=subs,ncol=voxels)
progress <- txtProgressBar(min = 0, max = voxels, style = 3)
for (i in 1:nvox){
	vox<-imat[,i]
	regfit<-lm(vox~age)
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

* Linear Regression
```
dm <-model.matrix(~age-1)
dm <-cbind(dm,1)
dm
   age  
1   69 1
2   75 1
3   75 1
4   78 1
5   82 1
6   84 1
7   71 1
8   72 1
9   58 1
10  66 1
11  86 1
12  81 1
13  52 1
14  90 1
15  80 1
16  92 1
17  71 1
18  73 1
19  75 1
20  70 1
21  73 1
22  61 1
23  61 1
24  62 1

conmat <-matrix(c(1,0),nrow=2,ncol=1)
```

* Multiple Regression
```
dm <-model.matrix(~age+mmse-1)
dm <-cbind(dm,1)
dm
   age mmse  
1   69   29 1
2   75   28 1
3   75   20 1
4   78   23 1
5   82   29 1
6   84   28 1
7   71   17 1
8   72   26 1
9   58   27 1
10  66   29 1
11  86   27 1
12  81   29 1
13  52   30 1
14  90   29 1
15  80   30 1
16  92   24 1
17  71   29 1
18  73   27 1
19  75   22 1
20  70   29 1
21  73   23 1
22  61   28 1
23  61   30 1
24  62   26 1

conmat <-matrix(0L,nrow=3,ncol=2)
conmat[,1] <-c(1,0,0) # association of age
conmat[,2] <-c(0,1,0) # association of MMSE
```

* ANOVA (one-way)
```
dm <-model.matrix(~sex-1)
dm <-cbind(dm,1)
dm
sexF sexM  
1     1    0 1
2     0    1 1
3     0    1 1
4     1    0 1
5     1    0 1
6     1    0 1
7     0    1 1
8     1    0 1
9     0    1 1
10    1    0 1
11    0    1 1
12    0    1 1
13    1    0 1
14    1    0 1
15    1    0 1
16    1    0 1
17    1    0 1
18    0    1 1
19    0    1 1
20    1    0 1
21    1    0 1
22    1    0 1
23    0    1 1
24    1    0 1

conmat <-matrix(0L,nrow=2,ncol=1)
conmat[,1] <-c(1,-1) # effect of sex
```

* ANOVA (two-way)
```
dm <-model.matrix(~sex:cdr-1)
dm <-cbind(dm,1)
dm
   sexF:cdrmild sexM:cdrmild sexF:cdrnotdem sexM:cdrnotdem sexF:cdrverymild
1             0            0              1              0                0
2             0            0              0              1                0
3             0            1              0              0                0
4             1            0              0              0                0
5             0            0              1              0                0
6             0            0              1              0                0
7             0            1              0              0                0
8             0            0              0              0                1
9             0            0              0              1                0
10            0            0              1              0                0
11            0            0              0              0                0
12            0            0              0              0                0
13            0            0              1              0                0
14            0            0              1              0                0
15            0            0              1              0                0
16            0            0              0              0                1
17            0            0              1              0                0
18            0            0              0              0                0
19            0            1              0              0                0
20            0            0              0              0                1
21            0            0              0              0                1
22            0            0              1              0                0
23            0            0              0              1                0
24            0            0              1              0                0
   sexM:cdrverymild  
1                 0 1
2                 0 1
3                 0 1
4                 0 1
5                 0 1
6                 0 1
7                 0 1
8                 0 1
9                 0 1
10                0 1
11                1 1
12                1 1
13                0 1
14                0 1
15                0 1
16                0 1
17                0 1
18                1 1
19                0 1
20                0 1
21                0 1
22                0 1
23                0 1
24                0 1

conmat <-matrix(0L,nrow=7,ncol=3)
conmat[,1] <-c(1, -1, 1, -1, 1, -1, 0) #Main effect of sex
conmat[,2] <-c(2, 2, -1, -1 , -1, -1, 0) #Clinical dementia rating using "not demented"
conmat[,3] <-c(2, -2, -1, 1, -1, 1, 0) # Interaction between two other contrasts
```

* ANCOVA
```
dm <-model.matrix(~age:sex-1)
dm <-cbind(dm,1)
dm
   age:sexF age:sexM  
1        69        0 1
2         0       75 1
3         0       75 1
4        78        0 1
5        82        0 1
6        84        0 1
7         0       71 1
8        72        0 1
9         0       58 1
10       66        0 1
11        0       86 1
12        0       81 1
13       52        0 1
14       90        0 1
15       80        0 1
16       92        0 1
17       71        0 1
18        0       73 1
19        0       75 1
20       70        0 1
21       73        0 1
22       61        0 1
23        0       61 1
24       62        0 1

conmat <-matrix(0L,nrow=3,ncol=3)
conmat[,1] <-c(1, 0, 0) #effect of age in females
conmat[,2] <-c(0, 1, 0) #effect of age in males
conmat[,3] <-c(1, -1, 0) #The ANCOVA design for comparing a covariate in groups
```

* Adding Controls
```
dm <-model.matrix(~sex+age-1)
dm <-cbind(dm,1)
dm
   sexF sexM age  
1     1    0  69 1
2     0    1  75 1
3     0    1  75 1
4     1    0  78 1
5     1    0  82 1
6     1    0  84 1
7     0    1  71 1
8     1    0  72 1
9     0    1  58 1
10    1    0  66 1
11    0    1  86 1
12    0    1  81 1
13    1    0  52 1
14    1    0  90 1
15    1    0  80 1
16    1    0  92 1
17    1    0  71 1
18    0    1  73 1
19    0    1  75 1
20    1    0  70 1
21    1    0  73 1
22    1    0  61 1
23    0    1  61 1
24    1    0  62 1

conmat <-matrix(0L,nrow=4,ncol=1)
conmat[,1] <-c(1, -1, 0, 0) #effect of sex when controlling for age
```

```
## Set up model
formula <-imat~var1-1
conmat <-matrix(c(1,0),nrow=2, ncol=1)
## fit model and extract important information
fit <-rft.lm(formula, conmat, mask, test="FALSE")
df <-fit$df
timg <-fit$tfields
fwhm <-fit$fwhm
## generate results
thresh <-rft.thresh(3, timg, .05, 150, fwhm, mask, df, "T", "voxel")
results <-rft.results(3, thresh[1], 100, fwhm, timg, mask, df, "T", thresh[2])
```



## Steps 3-6: Choosing a threshold and extracting important data

Fortunately, the final steps are intricately related and are consolidated into a handful of functions. An analysis can be completed using only `rft.thresh` and `rft.results` takes care of the rest. However, it can be beneficial to be aware of what functions are used within the `rft.results` functions.

** `rft.pcluster` - Calculates probability of obtaining a cluster of a certain size given a threshold value and total search area
** `ants.resels` - The resolutions per voxels of an image
** `ants.ec` - Determines the euler characteristic given the threshold to maxima for a cluster

```
thresh<-rft.thresh(D,timg,.05,150,fwhm,mask,rdf,"T","voxel")
results<-rft.results(D,thresh[1],100,fwhm,timg,mask,rdf,"T",thresh[2])
```


