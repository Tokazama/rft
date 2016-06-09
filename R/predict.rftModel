#' @param object
#' @param newImages
#' @param newdata data must have the same names as those specified in the previous formula
#' 
#'
#' @return 
#' \item{fit} {vector or matrix as above}
#' \item{se.fit} {standard error of predicted means}
#' \item{residual.scale} {residual standard deviations}
#' \item{dof} {degrees of freedom for residual}
#'
#'
#'
#' @export predict.rftModel
predict.rftModel <- function (object, newImages, newRH, newdata, contrast, conType,
                            isotropic = TRUE, k = 1, threshType = "pRFT", pval = .05,
                            pp = .001, n = 1, statdir, verbose = FALSE) {
  if (missing(newdata)) {
    data <- environment(RH)
  }
  call <- match.call()
  
  X <- model.matrix(RH, data = data)
  # standard error-------------------------------------------------------------
  r <- newImages - X %*% object$coefficients
  rMean <- colMeans(r)
  rSE <- sqrt(colSums((r - rMean)^2) / (n - 1)) / sqrt(n)
  totalSE <- sd(r) / sqrt(n)
  
  # re-estimate smoothness parameters------------------------------------------
  if (verbose)
    cat("Estimating FWHM. \n")
  smooth <- estSmooth(r, mask, rdf, scaleResid = FALSE, sample = sample, verbose = verbose)
  if (verbose)
    cat("Estimating resels. \n")
  v2r <- resels(mask, smooth$fwhm)
  
  connames <- names(contrastImagesTest$contrastImages)
  names(contrastsResultsTrain) <- connames
  names(contrastsResultsTest) <- connames
  
  structure(list(train = object, test = testModel), class = "predict.rftModel")
}

#' @export summary.predict.rftModel
summary.predict.rftModel <- function(object, contrastMatrix, statdir = NULL, controls = list(), verbose = FALSE) {
  # check arguments------------------------------------------------------------
  controlvals <- rftControl()
  if (!missing(control))
    controlvals[names(control)] <- control
  if (!missing(contrastMatrix))
    contrastMatrix <- object$contrastMatrix
  else {
    condim <- dim(contrastMatrix)
    if (is.null(condim))
      contrastMatrix <- matrix(contrastMatrix, 1)
  }
  ncon <- nrow(contrastMatrix)
  
  if (ncol(contrastMatrix) != object$volumes$np)
    stop("Each contrast length must be equal to the number of columns in the design matrix, including any intercepts. \n")
  if (is.null(rownames(contrastMatrix))) {
    conNames <- paste("contrastImage", 1:ncon, sep = "")
    rownames(contrastMatrix) <- conNames
  }
  else
    conNames <- rownames(contrastMatrix)
  if (is.null(rownames(contrastMatrix)))
    conNames <- paste("contrastImage", 1:ncon, sep = "")
  else
    conNames <- rownames(contrast)
  
  pred_summary <- matrix(0, ncon, 7)
  colnames(pred_summary) <- c("Train_SetLevel", "Test_SetLevel", "Train_Max_tstat", "Test_Max_tstat", "Train_Mean_SE", "Test_Mean_SE", "Dice")
  rownames(pred_summary) <- connames
  
  ans <- list()
  # compute contrasts----------------------------------------------------------
  for (i in 1:ncon) {
    c <- matrix(contrast[i, ], ncol = 1)
    # summarise each model independently
    trainsum <- summary(object$train, c, statdir = statdir, controls = controlvals, verbose = verbose)
    testsum  <- summary(object$test,  c, statdir = statdir, controls = controlvals, verbose = verbose)
    
    # compared results summary-------------------------------------------------
    pred_summary[i, 1] <- trainsum$set.level
    pred_summary[i, 2] <- testsum$set.level
    pred_summary[i, 3] <- max(testsum$contrastImage)
    pred_summary[i, 4] <- max(trainsum$contrastImage)
    pred_summary[i, 5] <- mean(voxelSE1[trainsum$clusterImage > 0])
    pred_summary[i, 6] <- mean(voxelSE2[tessumt$clusterImage > 0])
    pred_summary[i, 7] <- sum(trainsum$clusterImage > 0 &  testsum$clusterImage > 0) /
                          sum(trainsum$clusterImage > 0 |  testsum$clusterImage > 0)
    olimg <- makeImage(dim(object$train$dims$mask), 0)
    olimg[trainsum$clusterImage > 0] <- 1
    olimg[testsum$clusterImage > 0] <- 2
    olimg[trainsum$clusterImage > 0 && testsum$clusterImage > 0] <- 3
    
    
    ans <- lappend(ans, list(train = trainsum, test = testsum, overlapImage = olimg))
    names(ans)[i] <- conNames[i]
  }
  rownames(pred_summary) <- connames
  structure(list(contrasts = ans, comparison = pred_summary), class = "summary.predict.rftModel")
}

#' @export anova.predict.rftModel
anova.predict.rftModel <- function(object, contrastMatrix, statdir = NULL, controls = list(), verbose = FALSE) {
  # check arguments------------------------------------------------------------
  controlvals <- rftControl()
  if (!missing(control))
    controlvals[names(control)] <- control
  if (!missing(contrastMatrix))
    contrastMatrix <- object$contrastMatrix
  else {
    condim <- dim(contrastMatrix)
    if (is.null(condim))
      contrastMatrix <- matrix(contrastMatrix, 1)
  }
  ncon <- nrow(contrastMatrix)
  
  if (ncol(contrastMatrix) != object$volumes$np)
    stop("Each contrast length must be equal to the number of columns in the design matrix, including any intercepts. \n")
  if (is.null(rownames(contrastMatrix))) {
    conNames <- paste("contrastImage", 1:ncon, sep = "")
    rownames(contrastMatrix) <- conNames
  }
  else
    conNames <- rownames(contrastMatrix)
  if (is.null(rownames(contrastMatrix)))
    conNames <- paste("contrastImage", 1:ncon, sep = "")
  else
    conNames <- rownames(contrast)
  
  pred_summary <- matrix(0, ncon, 7)
  colnames(pred_summary) <- c("Train_SetLevel", "Test_SetLevel", "Train_Max_tstat", "Test_Max_tstat", "Train_Mean_SE", "Test_Mean_SE", "Dice")
  rownames(pred_summary) <- connames
  
  ans <- list()
  # compute contrasts----------------------------------------------------------
  for (i in 1:ncon) {
    c <- matrix(contrast[i, ], ncol = 1)
    # summarise each model independently
    trainsum <- anova(object$train, c, statdir = statdir, controls = controlvals, verbose = verbose)
    testsum  <- anova(object$test,  c, statdir = statdir, controls = controlvals, verbose = verbose)
    
    # compared results summary-------------------------------------------------
    pred_summary[i, 1] <- trainsum$set.level
    pred_summary[i, 2] <- testsum$set.level
    pred_summary[i, 3] <- max(testsum$contrastImage)
    pred_summary[i, 4] <- max(trainsum$contrastImage)
    pred_summary[i, 5] <- mean(voxelSE1[trainsum$clusterImage > 0])
    pred_summary[i, 6] <- mean(voxelSE2[tessumt$clusterImage > 0])
    pred_summary[i, 7] <- sum(trainsum$clusterImage > 0 &  testsum$clusterImage > 0) /
      sum(trainsum$clusterImage > 0 |  testsum$clusterImage > 0)
    olimg <- makeImage(dim(object$train$dims$mask), 0)
    olimg[trainsum$clusterImage > 0] <- 1
    olimg[testsum$clusterImage > 0] <- 2
    olimg[trainsum$clusterImage > 0 && testsum$clusterImage > 0] <- 3
    
    ans <- lappend(ans, list(train = trainsum, test = testsum, overlapImage = olimg))
    names(ans)[i] <- conNames[i]
  }
  rownames(pred_summary) <- connames
  structure(list(contrasts = ans, comparison = pred_summary), class = "anova.predict.rftModel")
}

#' @export print.summary.predict.rftModel
print.summary.predict.rftModel <- function(object) {
  print(round(as.data.frame(object$comparison), 2))
}

#' @export print.anova.predict.rftModel
print.anova.predict.rftModel <- function(object) {
  print(round(as.data.frame(object$comparison), 2))
}
