#' Color Gradient Bar
#'
#' @param lut color gradient function
#' @param min minimum color value to display (scale 1 - 100)
#' @param max maximum color value to display (scale 1 - 100)
#' @param side which side to place the gradient bar on. 
#'        Possible strings are:, 1 = left, 2 = right, 3 = bottom, 
#'        4 = top, 5 = center upright, 6 = center sideways
#' @param justExtent fraction indicating how far to justify from center
#' @param minLabel text to label minimum of gradient bar
#' @param maxLabel text to label maximum of gradient bar
#' 
#' 
gradientBar <- function(lut, min = 1, max = 100, side = 2, justExtent = .9, minLabel, maxLabel) {
  DIM <- dim(grid.layout())
  mycolors <- lut(100)
  if (side == 1 || side == 2 || side == 5) {
    mycolors <- rev(mycolors)
    if (side == 1) {
      grid::grid.raster(mycolors[min:max], x = DIM[1] * justExtent, y = .5, width = 0.075 * DIM[1], height = DIM[2] * 0.6)
      if (!missing(minLabel))
        grid::grid.text(minLabel, x = DIM[1] * justExtent, y = DIM[2] - (DIM[2] * .85))
      if (!missing(maxLabel))
        grid::grid.text(maxLabel, x = DIM[1] * justExtent, y = DIM[2] * 0.85)
    } else if (side == 2) {
      grid::grid.raster(mycolors[min:max], x = DIM[1] - DIM[1] * justExtent, y = .5, width = 0.075 * DIM[1], height = DIM[2] * 0.6)
      if (!missing(minLabel))
        grid::grid.text(minLabel, x = DIM[1] - DIM[1] * justExtent, y = DIM[2] - (DIM[2] * .85))
      if (!missing(maxLabel))
        grid::grid.text(maxLabel, x = DIM[1] - DIM[1] * justExtent, y = DIM[2] * 0.85)
    } else if (side == 5) {
      grid::grid.raster(mycolors[min:max], x = DIM[1] * .5, y = dim[2] * .5, width = 0.075 * DIM[1], height = DIM[2] * 0.6)
      if (!missing(minLabel))
        grid::grid.text(minLabel, x = DIM[1] * .5, y = DIM[2] - (DIM[2] * .85))
      if (!missing(maxLabel))
        grid::grid.text(maxLabel, x = DIM[1] * .5, y = DIM[2] * 0.85)
    }
  } else if (side == 3 || side == 4 || side == 6) {
    if (side == 3) {
      grid::grid.raster(t(mycolors[min:max]), x = 0.5, y = DIM[2] * justExtent, width = 0.6 * DIM[1], height = 0.075 * DIM[2])
      if (!missing(minLabel)) {
        minLabel_chars <- nchar(minLabel) * (DIM[1] * 0.01)
        grid::grid.text(minLabel, x = .2 * DIM[1] - minLabel_chars, y = DIM[2] * justExtent)
      }
      if (!missing(maxLabel))
        grid::grid.text(maxLabel, x = (.8 * DIM[1]) + DIM[1] * .02, y = DIM[2] * justExtent)
    } else if (side == 4) {
      grid::grid.raster(t(mycolors[min:max]), x = 0.5, y = DIM[2] - DIM[2] * justExtent, width = 0.6 * DIM[1], height = 0.075 * DIM[2])
      if (!missing(minLabel)) {
        minLabel_chars <- nchar(minLabel) * (DIM[1] * 0.01)
        grid::grid.text(minLabel, x = .2 * DIM[1] - minLabel_chars, y = DIM[2] - DIM[2] * justExtent)
      }
      if (!missing(maxLabel)) 
        grid::grid.text(maxLabel, x = (.8 * DIM[1]) + DIM[1] * .02, y = DIM[2] - DIM[2] * justExtent)
    } else if (side == 6) {
      grid::grid.raster(t(mycolors[min:max]), x = DIM[1] * 0.5, y = DIM[2] * .5, width = 0.6 * DIM[1], height = 0.075 * DIM[2])
      if (!missing(minLabel)) {
        minLabel_chars <- nchar(minLabel) * (DIM[1] * 0.01)
        grid::grid.text(minLabel, x = .2 * DIM[1] - minLabel_chars, y = DIM[2] - DIM[2] * justExtent)
      }
      if (!missing(maxLabel)) 
        grid::grid.text(maxLabel, x = DIM[1] * 5, y = DIM[2] - DIM[2] * justExtent)
    }
  }
}
