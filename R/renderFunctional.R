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
  img <- antsImageClone(funcimg)
  # Perona malik edge preserving smoothing
  # iMath(img, "PeronaMalike", iterations (ex. 10), conductance (ex. .5) 
  
  img <- smoothImage(funcimg, smoothval)
  mask <- antsImageClone(img)
  mask[mask != 0] <- 1
  mask <- iMath(mask, "FillHoles")
  emask <- iMath(mask, "ME", 1)
  surf <- as.array(mask - emask)
  xyz <- which(surf != 0, arr.ind = TRUE)
  func <- as.array(img)
  unique_vox <- length(unique(func))
  if (missing(max)) {
    if (unique_vox > 2^8)
      max <- 2^8
    else
      max <- unique_vox
  }
  func <- round((func - min(func)) / max(func - min(func)) * (max - min) + min)
  vox_val <- as.numeric(func, sur > 0)
  
  
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
  
  # render surface-------------------------------------------------------------
  if (verbose)
    cat("Computing surface \n")
  brain <- misc3d::contour3d(mask, level = surfval, alpha = alpha, draw = FALSE, smooth = 1, material = material, depth = depth)

  if (length(colorGradient > 1)) {
    n <- length(brain$v1[, 1])
    vox_vec <- rep(0, n)
    findColor <- function(x) {
      surf_img[floor(x[i, 1]), floor(x[i, 2]), floor(x[i, 2])]
    }
    brain$color <- apply(brain$v1, 1, findColor)
    brain$color2 <- apply(brain$v2, 1, findColor)
    brain$color.mesh <- apply(brain$v3, 1, findColor)
  } else {
    brain$color <- mycolors
    # surf_img$color2 <- mycolors[myv2]
    # surf_img$col.mesh <- mycolors[myv3]
  }
    
  brain$v1 <- antsTransformIndexToPhysicalPoint(img, brain$v1)
  brain$v2 <- antsTransformIndexToPhysicalPoint(img, brain$v2)
  brain$v3 <- antsTransformIndexToPhysicalPoint(img, brain$v3)
  
  .check3d()
  misc3d::drawScene.rgl(brain)
  brain
}

origin <- antsGetOrigin(img)
DIM <- dim(img)

as.antsImage(img[1:origin[1], 1:DIM[2], 1:DIM[3]])
coronal_post <- img[1:DIM[1], 1:origin[2], 1:DIM[3]]
coronal_post <- img[1:DIM[1], origin[2]:DIM[2], 1:DIM[3]]

left_sagittal <- img[1:origin[1], 1:DIM[2], 1:DIM[3]]
right_sagittal <- img[origin[1]:DIM[1], 1:DIM[2], 1:DIM[3]]

superior_axial 
inferior_axial





img <- smoothImage(funcimg, smoothval)
mask <- antsImageClone(img)
mask[mask != 0] <- 1
mask <- iMath(mask, "FillHoles")
emask <- iMath(mask, "ME", 1)
surf <- as.array(mask - emask)
xyz <- which(surf != 0, arr.ind = TRUE)
spheres3d(xyz[, 1], xyz[, 2], xyz[, 3])
