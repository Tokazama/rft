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
#' @param depth
#' @param verbose
#' mnit <- antsImageRead(getANTsRData('mni'))
#'
#'
#'
#'
#' @export renderFunctional
renderFunctional <- function(funcimg, surfval = .5, colorGradient = "rainbow", max, min = 0, upper = 1, lower = 1, alpha = 1, smoothval = 0, material = "metal", depth = .6, verbose = TRUE) {
  # mask <- antsImageClone(img)
  # if (max(img) != 1 | min(img) != 0)
  #   mask[mask !=0] <- 1
  # surf_mask <- as.array(mask - iMath(mask, "ME", 1))
  # surf_img <- surf_mask * as.array(img)
  # mycoord <- which(surf_img != 0, arr.ind = TRUE) # coordinate for each voxel
  
  img <- antsImageClone(funcimg)
  img <- smoothImage(funcimg, smoothval)
  surf_img <- as.array(img)
  mask <- surf_img
  mask[mask != 0] <- 1
  if (missing(max)) {
    if (unique_vox > 2^8)
      max <- 2^8
    else
      max <- unique_vox
  }
  surf_img <- round((surf_img - min(surf_img)) / max(surf_img - min(surf_img)) * (max - min) + min)
  
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
    mycolors <- colorGradient(max, alpha = 0)
  
  # white out upper or lower percentage of image
  if (lower < 1)
    mycolors[1:floor(lower * max)] <- "white"
  if (upper < 1)
    mycolors[ceiling(upper * max):max] <- "white"
  
  # render surface-------------------------------------------------------------
  if (verbose)
    cat("Computing surface \n")
  brain <- misc3d::contour3d(mask, level = surfval, alpha = alpha, draw = FALSE, smooth = 1, material = material, depth = depth)
  
  n <- length(brain$v1[,1])
  vox_vec <- rep(0, n)
  for (i in 1:n) {
    myv1[i] <- surf_img[brain$v1[i, 1], brain$v1[i, 2], brain$v1[i, 2]]
    myv2[i] <- surf_img[brain$v2[i, 1], brain$v2[i, 2], brain$v2[i, 2]]
    myv3[i] <- surf_img[brain$v3[i, 1], brain$v3[i, 2], brain$v3[i, 2]]
  }
  
  brain$color <- mycolors[myv1]
  # surf_img$color2 <- mycolors[myv2]
  # surf_img$col.mesh <- mycolors[myv3]
  
  brain$v1 <- antsTransformIndexToPhysicalPoint(img, brain$v1)
  brain$v2 <- antsTransformIndexToPhysicalPoint(img, brain$v2)
  brain$v3 <- antsTransformIndexToPhysicalPoint(img, brain$v3)
  
  .check3d()
  misc3d::drawScene.rgl(brain)
  brain
}
