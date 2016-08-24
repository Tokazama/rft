library(ANTsR)
library(h5)


home <- '/Volumes/SANDISK/datasets/ucsd/'
ucsd <- read.csv(paste(home, 'spreadsheets/ucsdWOna.csv', sep = ""))[, -1]

path2rft <- "Desktop/git/rft/R/"
source(paste(path2rft, "iClass.R", sep = ""))
source(paste(path2rft, "iModel.R", sep = ""))
source(paste(path2rft, "ilm.R", sep = ""))
source(paste(path2rft, "iUtils.R", sep = ""))
source(paste(path2rft, "iContrast.R", sep = ""))
source(paste(path2rft, "iFilter.R", sep = ""))
source(paste(path2rft, "rftResults.R", sep = ""))

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
rownames(contrastMatrix) <- c("pos", "neg")
fit1 <- summary(fit1, contrastMatrix, cthresh = 150)
report(fit1, paste(home, "report.pdf"))

# ANOVA----
fit2 <- ilm(y ~ Injury:Gender - 1)
contrastMatrix <- matrix(nrow = 5, ncol = 4)
                         # Male*OI  # Female*OI  # Male*TBI  # Female*TBI
contrastMatrix[1, ] <- c( 1,       -1,           1,         -1)   # the main effect of Gender
contrastMatrix[2, ] <- c(-1,        1,          -1,          1)
contrastMatrix[3, ] <- c( 1,        1,          -1,         -1)   # the main effect of Injury
contrastMatrix[4, ] <- c(-1,       -1,           1,          1)
contrastMatrix[5, ] <- c( 1,       -1,          -1,          1)   # the interaction between Gender and Injury

rownames(contrastMatrix) <- c("M > F", "M < F", "OI > TBI", "OI < TBI", "Gender x Injury")
fit2 <- anova(fit4, contrastMatrix, cthresh = 100)

# predict----

# predict.summary----

# predict.anova----
