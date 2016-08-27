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

mask <- getMask(antsImageRead(wblist[1]), cleanup = 0)
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
report(fit1, paste(home, "report.pdf"))

# ANOVA----
fit2 <- ilm(wb ~ Injury:Gender - 1, mydata)
contrastMatrix <- matrix(nrow = 5, ncol = 4)
                         # Female*OI  # Female*TBI  # Male*OI  # Female*TBI
contrastMatrix[1, ] <- c( 1,          1,           -1,         -1)   # the main effect of Gender
contrastMatrix[2, ] <- c(-1,         -1,            1,          1)
contrastMatrix[3, ] <- c( 1,         -1,            1,         -1)   # the main effect of Injury
contrastMatrix[4, ] <- c(-1,          1,           -1,          1)
contrastMatrix[5, ] <- c( 1,         -1,           -1,          1)   # the interaction between Gender and Injury

rownames(contrastMatrix) <- c("F > M", " F < M", "OI > TBI", "OI < TBI", "Gender x Injury")
fit2 <- anova(fit2, contrastMatrix, cthresh = 100, threshType = "cFDR")

# ANCOVA----
fit3 <- ilm(wb ~ WASI.IQ:Injury, mydat)

cm1 <- matrix(nrow = 2, ncol = 4)
              # (Intercept) # OI  # TBI
cm1[1, ] <- c(0,            1,    -1)
cm1[2, ] <- c(0,           -1,     1)
rownames(cm1) <- c("OI > TBI", "OI < TBI")
fit2 <- anova(fit3, cm1, threshType = "cFDR")

cm2 <- matrix(nrow = 4, ncol = 4)
              # (Intercept) # OI  # TBI
cm1[1, ] <- c(0,            1,    0)
cm1[2, ] <- c(0,           -1,    0)
cm1[3, ] <- c(0,            0,    1)
cm1[4, ] <- c(0,            0,   -1)
rownames(cm2) <- c("OI+", "OI-", "TBI+", "TBI-")
fit2 <- summary(fit3, cm2)

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
