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


# rftGlm
fm3 <- rftGlm(imat, ~AgeAtScan, mask, data = ptbp, intercept = TRUE, conmat,
              statdir = '/Users/zach8769/Desktop/', verbose = TRUE)

# summary
sumfm3 <- summary(fm3)

# anova
aovfm3 <- anova(fm3)

# predict

# predict.summary

# predict.anova
