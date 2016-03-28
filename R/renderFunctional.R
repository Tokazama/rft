#' Render Functional Surface Image
#' 
#' Creates functional surface images with color gradient options using misc3d.
#' 
#' @param ImageList image to render
#' @param surfval intensity level that defines isosurface
#' @param color may either be a color function (i.e. rainbow()) for each image or a single color for the entire image
#' @param max integer. upper end of image scale
#' @param upper fraction of upper values to use
#' @param lower fraction of lower values to use
#' @param alpha opacity
#' @param smoothval smoothing for image
#' @param material options are "dull", "shiny", "metal", or "default"
#' @param draw logical.If \code{TRUE} an rgl surface image is produced
#' @param depth 
#' @return a 3D surface image (if \code{draw = TRUE}) and list that can be passed to drawScene.rgl()
#' @example
#' /dontrun{
#' mnit <- antsImageRead(getANTsRData('mni'))
#' myscene <- brainView(list(mnit), color = "rainbow")
#' brainColors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), interpolate = c("spline"), space = "Lab")
#'}
#' @export renderFunctional
renderFunctional <- function(ImageList, color, max, upper, lower,
                             alpha, smoothval, material, depth = depth,
                             draw = TRUE) {
  # set/check parameters-------------------------------------------------------
  nimg <- length(ImageList)
  if (nimg > 1) {
    for (i in 2:nimg) {
      if (dim(ImageList[[i-1]]) != dim(ImageList[[i]]))
        stop("All images must have equal dimensions")
    }
  }
  if (missing(max)) {
    max <- rep(0, nimg)
  } else if (length(max) != nimg) {
    cat("Length of max doesn't equal number of images. Resetting max.")
    max <- rep(0, nimg)
  }
  if (missing(depth))
    depth <- rep(.6, nimg)
  if (missing(upper))
    upper <- rep(1, nimg)
  if (missing(lower))
    lower <- rep(1, nimg)
  if (missing(material))
    material <- rep("material", nimg)
  if (missing(alpha)) {
    if (nimg > 1)
      alpha <- rep(.8, nimg)
    else
      alpha <- 1
  }
  eps <- .Machine$double.eps
  DIM <- dim(ImageList[[1]])
  scene <- list()
  # big for loop begins--------------------------------------------------------
  for (ifunc in 1:nimg) {
    if (missing(smoothval))
      imgar <- as.array(ImageList[[ifunc]])
    else
      imgar <- as.array(smoothImage(ImageList[[ifunc]], smoothval))
    unique_vox <- length(unique(imgar))
    if (max[ifunc] == 0) {
      if (unique_vox > 100)
        max[ifunc] <- 100
      else
        max[ifunc] <- unique_vox
    }
    # enforce voxel range
    imgar <- round((imgar - min(imgar)) / max(imgar - min(imgar)) * (max[ifunc] - 1) + 1)
    # acquire colors
    if (class(color[ifunc]) == "function") {
      mycolors <- color(max[ifunc])
      # white out upper or lower percentage of image
      if (lower < 1)
        mycolors[1:floor(lower * max[ifunc])] <- "white"
      if (upper < 1)
        mycolors[ceiling(upper * max[ifunc]):max[ifunc]] <- "white"
      # render surface-------------------------------------------------------------
      progress <- txtProgressBar(min = 0, max = max[ifunc], style = 3)
      level_seq <- c(1:max[ifunc])
      # render each level of functional images
      for (i in 1:max[ifunc]) {
        if (any(imgar == level_seq[i])) {
          tmp <- array(0, dim = dim(imgar))
          tmp[imgar == level_seq[i]] <- 1
          brain <- misc3d::contour3d(tmp, level = 0, alpha = alpha[ifunc],
                                     smooth = FALSE, depth = depth,
                                     color = mycolors[i], draw = FALSE,
                                     material = material)
          brain$v1 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v1)
          brain$v2 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v2)
          brain$v3 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v3)
          scene <- lappend(scene, list(brain))
        }
        setTxtProgressBar(progress, i)
      }
      close(progress)
    } else {
      # for solid color structures
      imgar[imgar != 0] <- 1
      brain <- misc3d::contour3d(imgar, level = 0, alpha = alpha[ifunc],
                                 smooth = FALSE, depth = depth,
                                 color = mycolors, draw = FALSE,
                                 material = material)
      brain$v1 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v1)
      brain$v2 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v2)
      brain$v3 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v3)
      scene <- lappend(scene, list(brain))
    }
  }
  if (draw == "TRUE")
    misc3d::drawScene.rgl(scene)
  scene
}
