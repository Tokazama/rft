#' @param img image to render
#' @param colorGradient
#' \itemize{
#' \item{terrain: } {}
#' \item{heat: } {heat map}
#' \item{topo: } {topological}
#' \item{cm: } {}
#' \item{gs: } {grey scale}
#' \item{rainbow: } {}
#' \item{brain: } {}
#' }
#' @param customGradient
#' @param ncolor number of distinct colors to use (default = 2^8)
#' @param alpha opacity
#' mnit <- antsImageRead(getANTsRData('mni'))
#'
#'
#'
#'
#' @export renderSurface
renderSurface <- function(img, colorGradient = "rainbow", customGradient = NULL,
                          ncolor = 2^8, alpha = 1) {
  # isolate surface voxels-----------------------------------------------------
  mask <- antsImageClone(img)
  if (max(img) != 1 | min(img) != 0)
    mask[mask !=0] <- 1
  surf_mask <- as.array(mask - iMath(mask, "ME", 1))
  surf_img <- surf_mask * as.array(img)
  mycoord <- which(surf_img != 0, arr.ind = TRUE) # coordinate for each voxel
  mycoord2 <- antsTransformIndexToPhysicalPoint(img, mycoord)
  
  # linMap <- function(x, from, to) (x - min(x)) / max(x - min(x)) * (to - from) + from
  func_img <- round((surf_img - min(surf_img)) / max(surf_img - min(surf_img)) * ncolor)
  
  # create color lookup table--------------------------------------------------
  vox_unique <- unique(func_img)
  brain.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),
                                   interpolate = c("spline"), space = "Lab")
  if (missing(colorGradient) && vox_unique > 30)
    colorGradient = "heat"
  else if (colorGradient[1] == "terrain")
    mycolors <- terrain.colors(ncolor, alpha = 0)
  else if (colorGradient[1]  == "heat")
    mycolors <- heat.colors(ncolor, alpha = 0)
  else if (colorGradient[1]  == "topo")
    mycolors <- topo.colors(ncolor, alpha = 0)
  else if (colorGradient[1]  == "cm")
    mycolors <- cm.colors(ncolor, alpha = 0)
  else if (colorGradient[1] == "brain")
    mycolors <- brain.colors(ncolor, alpha = 0)
  else if (colorGradient[1] == "rainbow")
    mycolors <- rainbow(ncolor, alpha = 0)
  else if (colorGradient[1] == "gs")
    mycolors <- grey.colors(n, alpha = 0)
  else if (!is.null(customGradient))
    mycolors <- customGradient
  col <- mycolors[func_img]
  .check3d()
  rgl.surface(mycoord[, 1], mycoord[, 2], mycoord[, 3], color = col, alpha = alpha)
}

color_lut <- function(x, y, z) {
  mycolors[surf_img[x, y, z]]
}

