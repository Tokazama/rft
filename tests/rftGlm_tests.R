# read in data
filepath <- '/Volumes/SANDISK/datasets/ptbp/rft/warp/'
ifiles <- file.path(filepath, list.files(filepath)[1:33])
mask <- getMask(antsImageRead(ifiles[1]))
imat <- imagesToMatrix(ifiles, mask)
ptbp <- read.csv('/Volumes/SANDISK/datasets/ptbp/data/ptbp_summary_demographics.csv')
X <- model.matrix(~ptbp$AgeAtScan)
conmat <- matrix(c(0, 1, 0, -1), 2, 2, byrow = TRUE)

# rftModel
fm1 <- rftModel(imat, X, mask, contrastMatrix, diag(n), findVar = TRUE, findResels = FALSE, "rftGlm")


# update

# residuals

# rftREML

# rftGlm

# summary

# anova

# predict

# predict.summary

# predict.anova
