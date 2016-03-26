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
renderFunctional <- function(ilist, color, max, upper = 1, lower = 1, alpha = 1, slice = NULL, axis, smoothval, material = "metal") {
  DIM <- dim(ilist[[i]])
  nimg <- length(ilist)
  if (missing(max)) {
    max <- rep(0, nimg)
  } else if (length(max) != nimg) {
    cat("Length of max doesn't equal number of images. Resetting max.")
    max <- rep(0, nimg)
  }
  part1 <- list()
  part2 <- list()
  eps <- .Machine$double.eps
  for (i in 1:nimg) {
    if (is.null(slice)) {
      if (missing(smoothval))
        img <- as.array(ilist[[i]])
      else
        img <- smoothImage(ilist[[i]], smoothval)
      part1 <- lappend(part1, img)
    } else {
      if (missing(smoothval))
        img <- as.array(ilist[[i]])
      else
        img <- smoothImage(ilist[[i]], smoothval)
      if (axis == 1) {
        part1 <- lappend(part1, list(img[1:slice, 1:DIM[2], 1:DIM[3]]))
        part1[[i]][slice[1],,] <- 0
        part2 <- lappend(part1, list(img[slice:DIM[1], 1:DIM[2], 1:DIM[3]]))
        part2[[i]][,slice,] <- 0
      } else if (axis == 2) {
        part1 <- lappend(part1, list(img[1:DIM[1], 1:DIM[2], 1:slice]))
        part1[[i]][,, slice] <- 0
        part2 <- lappend(part2, list(img[1:DIM[1], 1:DIM[2], slice:DIM[3]]))
        part2[[i]][,, slice] <- 0
      } else if (axis == 3) {
        part1 <- lappend(part1, list(img[1:DIM[1], 1:slice, 1:DIM[3]]))
        part1[[i]][, slice,] <- 0
        part2 <- lappend(part2, list(img[1:DIM[1], slice:DIM, 1:DIM[3]]))
        part2[[i]][, slice,] <- 0
      }
    }
  }
  # Perona malik edge preserving smoothing
  # iMath(img, "PeronaMalike", iterations (ex. 10), conductance (ex. .5) 
  for (ifunc in 1:ilist) {
    unique_vox <- length(unique(ilist[[ifunc]]))
    if (missing(max[[ifunc]])) {
      if (unique_vox > 2^7)
        max[[ifunc]] <- 2^7
      else
        ax[[ifunc]] <- unique_vox
    }
    func <- round((ilist[[ifunc]] - min(ilist[[ifunc]])) / max(ilist[[ifunc]] - min(ilist[[ifunc]])) * (max[[ifunc]] - 1) + 1)

    # create color lookup table--------------------------------------------------
    if (verbose)
      cat("Creating color palette. \n")
    
    brain.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),
                                     interpolate = c("spline"), space = "Lab")
    if (colorGradient[1] == "terrain")
      mycolors <- terrain.colors(max[[ifunc]], alpha = 0)
    else if (colorGradient[1]  == "heat")
      mycolors <- heat.colors(max[[ifunc]], alpha = 0)
    else if (colorGradient[1]  == "topo")
      mycolors <- topo.colors(max[[ifunc]], alpha = 0)
    else if (colorGradient[1]  == "cm")
      mycolors <- cm.colors(max[[ifunc]], alpha = 0)
    else if (colorGradient[1] == "brain")
      mycolors <- brain.colors(max[[ifunc]], alpha = 0)
    else if (colorGradient[1] == "rainbow")
      mycolors <- rainbow(max[[ifunc]], alpha = 0)
    else if (colorGradient[1] == "gs")
      mycolors <- grey.colors(max[[ifunc]], alpha = 0)
    else
      mycolors <- colorGradient
    # white out upper or lower percentage of image
    if (lower < 1)
      mycolors[1:floor(lower * max[[ifunc]])] <- "white"
    if (upper < 1)
      mycolors[ceiling(upper * max[[ifunc]]):max[[ifunc]]] <- "white"

    # render surface-------------------------------------------------------------
    if (verbose)
      cat("Computing surface \n")
    if (verbose)
      progress <- txtProgressBar(min = 0, max = max[[ifunc]], style = 3)
    level_seq <- c(1:max[[ifunc]])
      for (i in 1:max[[ifunc]]) {
        if (any(ilist[[ifunc]] == level_seq[i])) {
          tmp <- array(0, dim = DIM)
          tmp[ ilist[[ifunc]] == level_seq[[i]] ] <- 1
          brain <- contour3d(tmp, level = 0, alpha = alpha[i], color = mycolors[i], draw = FALSE)
          brain$v1 <- antsTransformIndexToPhysicalPoint(img, brain$v1)
          brain$v2 <- antsTransformIndexToPhysicalPoint(img, brain$v2)
          brain$v3 <- antsTransformIndexToPhysicalPoint(img, brain$v3)
          scene <- c(scene, list(brain))
        }
      if (verbose)
        setTxtProgressBar(progress, i)
      }
    }
    if (verbose)
      close(progress)
    scene
  }
