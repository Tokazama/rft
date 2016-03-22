#' @param img image to render
#' @param surfval intensity level that defines isosurface
#' @param color
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
#' @param max integer; upper end of image scale
#' @param min integer; lower end of image scale
#' @param upper fraction of upper values to use
#' @param lower fraction of lower values to use
#' @param colorBelowUpper color for values below upper fraction
#' @param colorAboveLower color for values above lower fraction
#' @param alpha opacity
#' @param smoothval 
#' @param material
#' @param depth
#' @param add logical; if TRUE, add to current graph.
#' @param verbose
#' @return
#'
#'
#'
#' @export renderFunctional
renderFunctional <- function(funcimg, surfval = .5, color = "rainbow", max, min = 0, upper = 1, lower = 1, colorBelowUpper = "white", colorAboveLower = "white", alpha = 1, smoothval = 0, material = "metal", depth = .6, add = FALSE, verbose = TRUE) {
  # mask <- antsImageClone(img)
  # if (max(img) != 1 | min(img) != 0)
  #   mask[mask !=0] <- 1
  # surf_mask <- as.array(mask - iMath(mask, "ME", 1))
  # surf_img <- surf_mask * as.array(img)
  # set image parameters-------------------------------------------------------
  img <- antsImageClone(funcimg)
  img <- smoothImage(funcimg, smoothval)
  surf_img <- as.array(img)
  # default max parameters
  if (missing(max)) {
    unique_vox <- length(unique(surf_img))
    if (unique_vox > 2^8)
      max <- 2^8
    else
      max <- unique_vox
  }
  surf_img <- round((surf_img - min(surf_img)) / max(surf_img - min(surf_img)) * (max - min) + min)
  
  # create color palette--------------------------------------------------
  if (verbose)
    cat("Creating color palette. \n")
  if (colorGradient == "terrain" | colorGradient == "heat" | colorGradient == "topo" |
      colorGradient == "cm" | colorGradient == "brain" | colorGradient == "gs" | colorGradient == "rainbow") {
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
      mycolors <- colorGradient(max, alpha = 0)
      
    # white out upper or lower percentage of image
    if (lower < 1)
      mycolors[1:floor(lower * max)] <- colorAboveLower
    if (upper < 1)
      mycolors[ceiling(upper * max):max] <- colorBelowUpper
    findColor <- function(x, y, z) {
      voxval <- surf_img[x, y, z]
      mycolors[voxval]
    }
  } else {
    findColor <- color
  }
  
  # render surface-------------------------------------------------------------
  if (verbose)
    cat("Computing surface \n")
  mycoord <- which(surf_img > min, arr.ind = )
  de <- misc3d::kde3d(mycoord[,1], mycoord[,2], mycoord[,3], n = 40,
              lims = c(range(mycoord[,1]), range(mycoord[,2]), range(mycoord[,3])))
  brain <- misc3d::contour3d(de$d, .1, de$x, de$y, de$z, color = findColor,
                             alpha = alpha,  draw = FALSE, smooth = 1,
                             material = material, depth = depth)
  # use header to get correct orientation
  brain$v1 <- antsTransformIndexToPhysicalPoint(img, brain$v1)
  brain$v2 <- antsTransformIndexToPhysicalPoint(img, brain$v2)
  brain$v3 <- antsTransformIndexToPhysicalPoint(img, brain$v3)
  
  .check3d()
  misc3d::drawScene.rgl(brain, add = add)
  brain
}
