#' ICSNP has hotelling's T-test
#' @export rftGlm

rftGlm <- function(formula, mask, data, contrast, conType, statdir = NULL, sample = NULL, verbose = TRUE,...) {
  if (missing(data))
    data <- environment(formula)
  cl <- match.call()
  X <- model.matrix(formula, data = data)
  y <- model.response(model.frame(formula, data = data))

  v <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  edf <- p - 1
  rdf <- n - p
  
  b   <- matrix(0, nrow = p, ncol = v) # coefficients
  rownames(b) <- colnames(X)
  rss <- matrix(0, nrow = 1, ncol = v) # residual sum of squares
  r   <- matrix(0, nrow = n, ncol = v) # residuals
  
  # fit each voxel-------------------------------------------------------------
  if (verbose) {
    cat("Fitting model. \n")
    progress <- txtProgressBar(min = 0, max = v, style = 3)
  }
  for (i in 1:v) {
    fit <- glm.fit(X, y[, i], ...)
    b[, i]   <- fit$coefficients
    rss[, i] <- sum(fit$residuals^2)
    r[, i]   <- fit$residuals / sqrt(rss[, i] / rdf) # scale residuals with mean residual sum of squares
    if (verbose)
      setTxtProgressBar(progress, i)
  }
  if (verbose)
    close(progress)

  # estimate fwhm/resels-------------------------------------------------------
  if (!is.null(statdir)) {
    rl <- paste(statdir, "residual_", 1:n, ".nii.gz", sep = "")
    if (verbose)
      cat("Writing residual images to director. \n")
    for (i in 1:n) {
      antsImageWrite(makeImage(mask, r[i, ]), rl[i])
    }
    r <- rl
  }
  if (verbose)
    cat("Estimating FWHM. \n")
  smooth <- estSmooth(r, mask, rdf, scaleResid = FALSE, sample = sample, verbose = verbose)
  if (verbose)
    cat("Estimating resels. \n")
  v2r <- resels(mask, smooth$fwhm)
  
  if (missing(contrast))
    contrast <- NULL
  if (missing(conType))
    conType <- NULL
  # create rft object-----------------------------------------------------
  z <- structure(list(call = cl,
                      coefficients = b,
                      residuals = r,
                      rss = rss,
                      RPVImage = smooth$RPVImage,
                      fwhm = smooth$fwhm,
                      resels = v2r,
                      mask = antsImageClone(mask),
                      df = c(edf, rdf),
                      dm = X,
                      contrast = contrast,
                      contrast.type = conType,
                      voxels = v,
                      statdir = statdir,
                      method = "rftGlm"),
                 class = "rft")
  z
}
