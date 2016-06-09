# read in data
filepath <- '/Volumes/SANDISK/datasets/ptbp/rft/warp/'
ifiles <- file.path(filepath, list.files(filepath)[1:33])
mask <- getMask(antsImageRead(ifiles[1]))
imat <- imagesToMatrix(ifiles, mask)
ptbp <- read.csv('/Volumes/SANDISK/datasets/ptbp/data/ptbp_summary_demographics.csv')
X <- model.matrix(~ptbp$AgeAtScan)
conmat <- matrix(c(0, 1, 0, -1), 2, 2, byrow = TRUE)

# rftModel
fm1 <- rftModel(X, imat, mask, conmat)

# rftREML
remlparams <- rftREML(rftmod$varParams$Cy, X, rftmod$varParams$V)

# update
fm2 <- update(fm1, V = remlparams$V)

# residuals


# rftGlm

# summary

# anova

# predict

# predict.summary

# predict.anova
