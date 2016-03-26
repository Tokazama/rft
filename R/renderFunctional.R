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
#' @param max integer. upper end of image scale
#' @param upper fraction of upper values to use
#' @param lower fraction of lower values to use
#' @param alpha opacity
#' @param smoothval
#' mnit <- antsImageRead(getANTsRData('mni'))
#'
#'
#'
#'
#' @export renderFunctional
renderFunctional <- function(ImageList, color, max, upper = 1, lower = 1,
                             alpha = 1, slice = NULL, axis, smoothval,
                             material = "metal") {
  # set/check parameters-------------------------------------------------------
  nimg <- length(ImageList)
  if (nimg > 1) {
    for (i in 2:nimg) {
      if (dim(ImageList[[i-1]]) != dim(ImageList[[i]]))
        stop("All images must have equal dimensions")
    }
  }
  if (missing(max)) {
    if (is.null(slice))
      max <- rep(0, nimg)
    else
      max <- rep(0, 2*nimg)
  } else if (length(max) != nimg) {
    cat("Length of max doesn't equal number of images. Resetting max.")
    if (is.null(slice))
      max <- rep(0, nimg)
    else
      max <- rep(0, 2*nimg)
  }
  ilist <- list()
  eps <- .Machine$double.eps
  DIM <- dim(ImageList[[1]])
  # slice/array images---------------------------------------------------------
  for (i in 1:nimg) {
    if (is.null(slice)) {
      if (missing(smoothval))
        ilist <- lappend(ilist, as.array(ImageList[[i]]))
      else
        # Perona malik edge preserving smoothing
        # iMath(img, "PeronaMalike", iterations (ex. 10), conductance (ex. .5) 
        ilist <- lappend(ilist, as.array(smoothImage(ImageList[[i]], smoothval)))
    } else {
      if (missing(smoothval))
        img <- as.array(ImageList[[i]])
      else
        img <- smoothImage(ImageList[[i]], smoothval)
      if (axis == 1) {
        ilist <- lappend(part1, list(img[1:slice, 1:DIM[2], 1:DIM[3]]))
        ilist[[i]][slice[1],,] <- 0
        ilist <- lappend(part1, list(img[slice:DIM[1], 1:DIM[2], 1:DIM[3]]))
        ilist[[i]][,slice,] <- 0
      } else if (axis == 2) {
        ilist <- lappend(part1, list(img[1:DIM[1], 1:DIM[2], 1:slice]))
        ilist[[i]][,, slice] <- 0
        ilist <- lappend(part2, list(img[1:DIM[1], 1:DIM[2], slice:DIM[3]]))
        ilist[[i]][,, slice] <- 0
      } else if (axis == 3) {
        ilist <- lappend(part1, list(img[1:DIM[1], 1:slice, 1:DIM[3]]))
        ilist[[i]][, slice,] <- 0
        ilist <- lappend(part2, list(img[1:DIM[1], slice:DIM, 1:DIM[3]]))
        ilist[[i]][, slice,] <- 0
      }
    }
  }
  scene <- list()
  nimg <- length(ilist)
  # big for loop begins--------------------------------------------------------
  for (ifunc in 1:nimg) {
    unique_vox <- length(unique(ilist[[ifunc]]))
    if (max[ifunc] == 0) {
      if (unique_vox > 2^7)
        max[ifunc] <- 2^7
      else
        max[ifunc] <- rainbow(nimg)[ifunc]
    }
    # enforce voxel range
    func <- round((ilist[[ifunc]] - min(ilist[[ifunc]])) / max(ilist[[ifunc]] - min(ilist[[ifunc]])) * (max[ifunc] - 1) + 1)

    # acquire colors
    if (color[ifunc] == "terrain" | color[ifunc] == "heat" | color[ifunc] == "topo" |
        color[ifunc] == "cm" | color[ifunc] == "brain" | color[ifunc] == "rainbow" |
        color[ifunc] == "gs" | class(color[ifunc]) == "function") {
      brain.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),
                                       interpolate = c("spline"), space = "Lab")
      if (color[ifunc] == "terrain")
        mycolors <- terrain.colors(max[ifunc], alpha = 0)
      else if (color[ifunc]  == "heat")
        mycolors <- heat.colors(max[ifunc], alpha = 0)
      else if (color[ifunc]  == "topo")
        mycolors <- topo.colors(max[ifunc], alpha = 0)
      else if (color[ifunc]  == "cm")
        mycolors <- cm.colors(max[ifunc], alpha = 0)
      else if (color[ifunc] == "brain")
        mycolors <- brain.colors(max[ifunc], alpha = 0)
      else if (color[ifunc] == "rainbow")
        mycolors <- rainbow(max[ifunc], alpha = 0)
      else if (color[ifunc] == "gs")
        mycolors <- grey.colors(max[ifunc], alpha = 0)
      else if (class(color[ifunc]) == "function")
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
          if (any(func == level_seq[i])) {
            tmp <- array(0, dim = dim(func))
            tmp[func == level_seq[i]] <- 1
            brain <- misc3d::contour3d(tmp, level = 0, alpha = alpha[ifunc], color = mycolors[i], draw = FALSE, material = material)
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
      tmp <- array(0, dim = dim(func))
      tmp[func != 0] <- 1
      brain <- misc3d::contour3d(tmp, level = 0, alpha = alpha[ifunc], color = mycolors, draw = FALSE, material = material)
      brain$v1 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v1)
      brain$v2 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v2)
      brain$v3 <- antsTransformIndexToPhysicalPoint(ImageList[[ifunc]], brain$v3)
      scene <- lappend(scene, list(brain))
    }
  }
  misc3d::drawScene.rgl(scene)
  scene
  }
