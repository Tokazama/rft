library(ANTsR)
library(h5)


home <- '/Volumes/SANDISK/datasets/ucsd/'
ucsd <- read.csv(paste(home, 'spreadsheets/ucsdWOna.csv', sep = ""))[, -1]

path2rft <- "Desktop/git/rft/R/"
source(paste(path2rft, "iClass.R", sep = ""))
source(paste(path2rft, "iModel.R", sep = ""))
source(paste(path2rft, "ilm.R", sep = ""))
source(paste(path2rft, "iREML.R", sep = ""))
source(paste(path2rft, "iUtils.R", sep = ""))
source(paste(path2rft, "rftResults.R", sep = ""))
source(paste(path2rft, "statFieldThresh.R", sep = ""))


# Read in VBM----
wblist <- c()
boolwb <- rep(FALSE, nrow(ucsd))
for (i in 1:nrow(ucsd)) {
  tmppath <- paste(home, "warp/", ucsd$file_names[i], ".nii.gz", sep = "")
  if (file.exists(tmppath)) {
    wblist <- c(wblist, tmppath)
    boolwb[i] <- TRUE
  }
}

mask <- antsImageRead("/Volumes/SANDISK/datasets/ucsd/template/T_template0_BrainCerebellumExtractionMask.nii.gz")
imat <- imagesToMatrix(wblist, mask)
wb <- iGroup(imat, "wb", mask, modality = "VBM", checkMask = TRUE)
iGroupWrite(wb, paste(home, "wb.h5", sep = ""))

# Read in CT----
ctlist <- c()
boolct <- rep(FALSE, nrow(ucsd))
for (i in 1:nrow(ucsd)) {
  tmppath <- paste(home, "/ct/", ucsd$file_names[i], ".nii.gz", sep = "")
  if (file.exists(tmppath)) {
    ctlist <- c(ctlist, tmppath)
    boolct[i] <- TRUE
  }
}

mask <- getMask(antsImageRead(ctlist[1]), cleanup = 0)
ct <- imagesToMatrix(ctlist, mask) %>% iGroup("ct", mask, modality = "CT", filename = paste(home, "ct.h5", sep = ""))

# Create iData----
ct <- iGroupRead(paste(home, "ct.h5", sep = ""))
wb <- iGroupRead(paste(home, "wb.h5", sep = ""))


(mydata <- iData(list(ct, wb), list(boolct, boolwb), ucsd))
iDataWrite(mydata, paste(home, "iData", sep = ""))
(mydata <- substract(mydata, "wb"))
mydata <- add(mydata, wb, boolwb)

# Regression----
mydata <- iDataRead(paste(home, "iData", sep = ""))
fit1 <- ilm(wb ~ Age, mydata)
                           # Intercept # Age
contrastMatrix <- matrix(c(0,          1,    # positive correlation
                           0,         -1),   # negative correlation
                         2, 2, byrow = TRUE)
rownames(contrastMatrix) <- c("Age +", "Age -")
fit1 <- summary(fit1, contrastMatrix, cthresh = 150)
report(fit1, paste(home, "Age", sep = ""))

# ANOVA----
fit2 <- ilm(wb ~ Injury:Gender - 1, mydata)
contrastMatrix <- matrix(nrow = 5, ncol = 4)
                         # Female*OI  # Female*TBI  # Male*OI  # Female*TBI
contrastMatrix[1, ] <- c( 1,          1,           -1,         -1)    # the main effect of Gender
contrastMatrix[2, ] <- c(-1,         -1,            1,          1)
contrastMatrix[3, ] <- c( 1,         -1,            1,         -1)   # the main effect of Injury
contrastMatrix[4, ] <- c(-1,          1,           -1,          1)
contrastMatrix[5, ] <- c( 1,         -1,           -1,          1)   # the interaction between Gender and Injury

rownames(contrastMatrix) <- c("F > M", " F < M", "OI > TBI", "OI < TBI", "Gender x Injury")
fit2 <- anova(fit2, contrastMatrix, cthresh = 100, threshType = "cFDR")

# ANCOVA----
fit3 <- ilm(wb ~ WASI.IQ:Injury, mydata)

cm1 <- matrix(nrow = 2, ncol = 3)
              # (Intercept) # OI  # TBI
cm1[1, ] <- c(0,            1,    -1)
cm1[2, ] <- c(0,           -1,     1)
rownames(cm1) <- c("OI > TBI", "OI < TBI")
fit2 <- anova(fit3, cm1, threshType = "cFDR")

cm2 <- matrix(nrow = 4, ncol = 3)
              # (Intercept) # OI  # TBI
cm2[1, ] <- c(0,            1,    0)
cm2[2, ] <- c(0,           -1,    0)
cm2[3, ] <- c(0,            0,    1)
cm2[4, ] <- c(0,            0,   -1)
rownames(cm2) <- c("OI+", "OI-", "TBI+", "TBI-")
fit2 <- summary(fit3, cm2)

report(fit2)

# REML----
fit_reml <- ilm(wb ~ WASI.IQ, mydata, optim = "REML")

fit_iwls <- ilm(wb ~ WASI.IQ, mydata, optim = "IWLS")

cm3 <- matrix(nrow = 2, ncol = 2)
cm3[1, ] <- c(0,  1)
cm3[2, ] <- c(0, -1)

fit_reml <- summary(fit_reml, cm3)

fit_iwls <- summary(fit_iwls, cm3)






# plot----
# iData plotting
plot(mydata, "wb", "Age", fit1@C[[1]]$clusterImage)  # should plot specific cluster
plot(mydata)  # should just plot demographics

# iModel plotting
plot(fit1, "Age +", 1)  # should plot specific cluster

# predict----

# predict.summary----

# predict.anova----


# test a bunch of regression variables
vars <- ucsd[2:24]
varbool <- c()
for (i in 1:length(colnames(vars)))
  varbool <- c(varbool, is.factor(vars[, i]))

vars <- colnames(vars[!varbool])

progress <- txtProgressBar(min = 0, max = length(vars), style = 3)
for (i in seq_len(length(vars))) {
  
  ff <- as.formula(paste("wb ~ ", vars[i], sep = ""))
  fit <- ilm(ff, mydata)
  # Intercept # Variable
  contrastMatrix <- matrix(c(0,          1,    # positive correlation
                             0,         -1),   # negative correlation
                           2, 2, byrow = TRUE)
  rownames(contrastMatrix) <- c("Correlation +", "Correlation -")
  fit <- summary(fit, contrastMatrix)
  if (class(fit@C[[1]]$results) != "character" && class(fit@C[[2]]$results) != "character")
    report(fit, paste(home, vars[i], sep = ""))
  setTxtProgressBar(progress, i)
}
close(progress)

# perform a bunch of ANCOVA

for (i in seq_len(length(vars))) {
  
  ff <- as.formula(paste("wb ~ Injury:", vars[i], sep = ""))
  fit <- ilm(ff, mydata)

  cm1 <- matrix(nrow = 2, ncol = 3)
  # (Intercept) # OI  # TBI
  cm1[1, ] <- c(0,            1,    -1)
  cm1[2, ] <- c(0,           -1,     1)
  rownames(cm1) <- c("OI > TBI", "OI < TBI")
  fit <- anova(fit, cm1, threshType = "cFDR")
  
  cm2 <- matrix(nrow = 4, ncol = 3)
  # (Intercept) # OI  # TBI
  cm2[1, ] <- c(0,            1,    0)
  cm2[2, ] <- c(0,           -1,    0)
  cm2[3, ] <- c(0,            0,    1)
  cm2[4, ] <- c(0,            0,   -1)
  rownames(cm2) <- c("OI+", "OI-", "TBI+", "TBI-")
  fit <- summary(fit, cm2)
  report(fit, paste(home, "reports/", vars[i], sep = ""))
}

dir.create(paste(home, "reports/Sex", sep = ""))
for (i in seq_len(length(vars))) {
  # control for sex
  cat(vars[i], "\n\n")
  
  ff <- as.formula(paste(wb ~ Injury:Gender:, vars[i], sep = ))
  fit <- ilm(ff, mydata)
  
  cmf <- matrix(0, 9, 5)
  
  rownames(cmf) <- c("OI > TBI", "OI < TBI", "Female > Male", "Female < Male", "OI.Female > TBI.Female", "OI.Female < TBI.Female", "OI.Male > TBI.Male", "OI.Male < TBI.Male", "Injury.X.Sex")
  
  # Intercept, Female.OI, Female.TBI, Male.OI, Male.TBI
  
  cmf[1, ] <- c(0,  1, -1,  1, -1)
  cmf[2, ] <- c(0, -1,  1, -1,  1)
  cmf[3, ] <- c(0,  1,  1, -1, -1)
  cmf[4, ] <- c(0, -1, -1,  1,  1)
  cmf[5, ] <- c(0,  1, -1,  0,  0)
  cmf[6, ] <- c(0, -1,  1,  0,  0)
  cmf[7, ] <- c(0,  0,  0,  1, -1)
  cmf[8, ] <- c(0,  0,  0, -1,  1)
  cmf[9, ] <- c(0,  1, -1, -1,  1)
  fit <- anova(fit, cmf)
  
  cmt <- matrix(0,16, 5)
  rownames(cmt) <- c("Female.OI+", "Female.OI-", "Female.TBI+", "Female.TBI-", "Male.OI+", "Male.OI-", "Male.TBI+", "Male.TBI-", "Female+", "Female-", "Male+", "Male-", "OI+", "OI-", "TBI+", "TBI-")
  cmt[1, ] <- c(0,  1, 0, 0, 0)
  cmt[2, ] <- c(0, -1, 0, 0, 0)
  cmt[3, ] <- c(0, 0,  1, 0, 0)
  cmt[4, ] <- c(0, 0, -1, 0, 0)
  cmt[5, ] <- c(0, 0, 0,  1, 0)
  cmt[6, ] <- c(0, 0, 0, -1, 0)
  cmt[7, ] <- c(0, 0, 0, 0,  1)
  cmt[8, ] <- c(0, 0, 0, 0, -1)
  cmt[9, ]   <- c(0,  .5,  .5, 0, 0)
  cmt[10, ] <- c(0, -.5, -.5, 0, 0)
  cmt[11, ] <- c(0, 0, 0,  .5,   .5)
  cmt[12, ] <- c(0, 0, 0, -.5, -.5)
  cmt[13, ] <- c(0, .5, 0, .5, 0)
  cmt[14, ] <- c(0, -.5, 0, -.5, 0)
  cmt[15, ] <- c(0, 0, .5, 0, .5)
  cmt[16, ] <- c(0, 0, -.5, 0, -.5)
  fit <- summary(fit, cmt)
  
  report(fit, paste(home, "reports/Sex/", vars[i], sep = ""))
}


 




