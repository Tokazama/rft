#' @param img image to render
#' @param surfval intensity level that defines isosurface
#' @param colorGradient
#' \itemize{
#' \item{terrain: } {}
#' \item{heat: } {heat map}
#' \item{topo: } {topological}
#' \item{cm: } {}
#' \item{gs: } {grey scale}
#' \item{rainbow: } {}
#' \item{brain: } {}
#' \item{colorRamp: } {custom colorRamp function}
#' }
#' @param max integer. upper end of image scale
#' @param min integer. lower end of image scale
#' @param upper fraction of upper values to use
#' @param lower fraction of lower values to use
#' @param alpha opacity
#' @param smoothval
#' @param material
#' mnit <- antsImageRead(getANTsRData('mni'))
#'
#'
#'
#'
#' @export renderFunctional
renderFunctional <- function(funcimg, colorGradient = "rainbow", max, min = 1, upper = 1, lower = 1, alpha = 1, smoothval, material = "metal") {
  if (missing(smoothval))
    img <- antsImageClone(funcimg)
  else
    img <- smoothImage(funcimg, smoothval)

  # Perona malik edge preserving smoothing
  # iMath(img, "PeronaMalike", iterations (ex. 10), conductance (ex. .5) 
  func <- as.array(img)
  unique_vox <- length(unique(func))
  if (missing(max)) {
    if (unique_vox > 2^8)
      max <- 2^8
    else
      max <- unique_vox
  }
  func <- round((func - min(func)) / max(func - min(func)) * (max - min) + min)
  vox_vec <- as.numeric(func, func > 0)
  
  # create color lookup table--------------------------------------------------
  if (verbose)
    cat("Creating color palette. \n")
  
  brain.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),
                                   interpolate = c("spline"), space = "Lab")
  if (colorGradient[1] == "terrain") 
    mycolors <- terrain.colors(max, alpha = 0)
  else if (colorGradient[1]  == "heat") 
    mycolors <- heat.colors(max, alpha = 0)
  else if (colorGradient[1]  == "topo") 
    mycolors <- topo.colors(max, alpha = 0)
  else if (colorGradient[1]  == "cm") 
    mycolors <- cm.colors(max, alpha = 0)
  else if (colorGradient[1] == "brain") 
    mycolors <- brain.colors(max, alpha = 0)
  else if (colorGradient[1] == "rainbow") 
    mycolors <- rainbow(max, alpha = 0)
  else if (colorGradient[1] == "gs") 
    mycolors <- grey.colors(max, alpha = 0)
  else 
    mycolors <- colorGradient
  # white out upper or lower percentage of image
  if (lower < 1)
    mycolors[1:floor(lower * max)] <- "white"
  if (upper < 1)
    mycolors[ceiling(upper * max):max] <- "white"
  cols <- mycolors[vox_vec]
  
  # render surface-------------------------------------------------------------
  if (verbose)
    cat("Computing surface \n")
  if (verbose)
        progress <- txtProgressBar(min = 0, max = max, style = 3)
  level_seq <- c(min:max)
  for (i in 1:max) {
    eps <- .Machine$double.eps
    tmp <- array(FALSE, dim = (dim(img)))
    mlev <- level_seq
    if (i == max)
      nlev <- max + eps
    else
      nlev <- level_seq[i + 1]
    tmp[ func >= mlev & func < nlev ] <- 1
    if (sum(tmp != 0, na.rm=TRUE) == 0) {
      warning("No contour to make")
      next
    } else {
    brain <- contour3d(tmp, level = 0, alpha = alpha[i], color = cols[i], draw = FALSE)
    brain$v1 <- antsTransformIndexToPhysicalPoint(img, brain$v1)
    brain$v2 <- antsTransformIndexToPhysicalPoint(img, brain$v2)
    brain$v3 <- antsTransformIndexToPhysicalPoint(img, brain$v3)
    }
    scene <- c(scene, list(brain))
    if (verbose)
      setTxtProgressBar(progress, i)
  }
  if (verbose)
      close(progress)
  scene
}
