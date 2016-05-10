#' @param x  statistical field in matix form (1 x nvox)
#'
#' @references
#' Barnes et al., (2013) Set-level threshold-free tests on the intrinsic volumes of SPMs
#'
#' @export rftSetlevel

rftSetLevel <- function(object, mask, k, interval = .5) {
  p <- dim(coefficients(object))[1]
  r <- residuals(object)
  nsub <- dim(r)[1]
  edf <- p - 1
  rdf <- nsub - p
  myfwhm <- estSmooth(residuals(object), mask, rdf)
  myresel <- resels(mask, myfwhm$fwhm)
  # standardize residuals and create test image--------------------------------
  r <- r / (colSums(r) / rdf)
  r <- rbind(colMeans(r), r)
  
  # set gradient and loop constants--------------------------------------------
  glevels <- seq(min(r), max(r), by = interval)
  vox2res <- 1 / prod(r)
  k <- k * vox2res
  iter <- length(glevels)
  nimg <- dim(r)[1]
  LKC <- matrix(nrow = nimg, ncol = iter)
  
  # get set levels-------------------------------------------------------------
  for (i in 1:iter) {
    for (j in 1:nimg) {
      if (setLevels[i, 1] > 0)
        clust <- labelClusters(makeImage(mask, object$residuals[j,]), k, glevels[i], Inf)
      else if (glevels[i] < 0)
        clust <- labelClusters(makeImage(mask, object$residuals[j,]), k, -Inf, glevels[i])
      nclus <- length(unique(clust[clust > 0]))
      LKC[j, i] <- rftPval(mask@dimension, nclus, k, glevels[i], n, myresel, c(dfe, rdf), fieldType = "Z")$Pcor
    }
  }
  
  # multivariate analysis------------------------------------------------------
  X <- rbind(c(1, 0), cbind(rep(0, n), rep(1, n)))
  fit <- anova(lm(LKC ~ X))
  sumfit <- summary(fit, test = "Wilks")
  
  # empericial mean
  meanECtest <- mean(allLKCreg) * allpju_test
  # EC sd of original images
  sdECtest <- sd(allLKCreg) * allpju_test
  # EC profile based on image smoothness
  meanECtest_resel <- LKCresel * allpju_test
  
  # probably should just use generic plot.rftSetLevel for this
  # plot results---------------------------------------------------------------
  
  # title = Probability that this is a random field  p<%3.4f ',CVA.p
  # legend
  # xlabel = threshold
  # ylabel = EC
  # observed EC for s field', test_stat
  # random filed based on regression, test_stat
  #random field based on smoothness, test_stat
  }
