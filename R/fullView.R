#'
#' @param axis the axis to slice (1 , 2 or 3)
#' @param slice 
#' bm <-getMask( antsImageRead( getANTsRData("ch2") ) )
#' brain <-renderSurfaceFunction( surfimg =list( bm ) , alphasurf=0.1, smoothsval = 1.5 )
#' fullView(tempfile( fileext='.png'))
#'
#' @export fullView
fullView <- function(SurfaceImage, FunctionalImage, views, slice, axis, statdir){
  # create indicies of rotation--------------------------
  x90p <- rotationMatrix(pi/2, 1, 0, 0)
  x90n <- rotationMatrix(-pi/2, 1, 0, 0)
  x180p <- rotationMatrix(pi, 1, 0, 0)
  x180n <- rotationMatrix(pi, 1, 0, 0)
  
  y90p <- rotationMatrix(pi/2, 0, 1, 0)
  y90n <- rotationMatrix(-pi/2, 0, 1, 0)
  y180p <- rotationMatrix(pi, 0, 1, 0)
  y180n <- rotationMatrix(-pi, 0, 1, 0)
  
  z90p <- rotationMatrix(pi/2, 0, 0, 1)
  z90n <- rotationMatrix(-pi/2, 0, 0, 1)
  z180p <- rotationMatrix(pi, 0, 0, 1)
  z180n <- rotationMatrix(-pi, 0, 0, 1)
  
  # extract each viewpoint
  if (any(views == "anterior")) {
    rgl::par3d(userMatrix = x90n) # anterior view
    rgl::rgl.snapshot(paste(statdir, "anterior.png", sep = ""), fmt = "png", top = TRUE)
    aa <- png::readPNG(paste(statdir, "anterior.png", sep = ""))
    png(paste(statdir, "anterior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
    grid::grid.raster(aa)
    grid::grid.text("Anterior", .5, .95)
    grid::grid.text("L", .9, .5) 
    grid::grid.text("R", .1, .5)
    grid::grid.text("S", .5, .9)
    grid::grid.text("I", .5, .1)
    dev.off()
  }
  if (any(views == "lef")) {
    rgl::par3d(userMatrix = x90n %*% z90n) # left view
    rgl::rgl.snapshot(paste(statdir, "left.png", sep = ""), fmt = "png", top = TRUE)
    aa <- png::readPNG(paste(statdir, "left.png", sep = ""))
    png(paste(statdir, "left.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
    grid::grid.raster(aa)
    grid::grid.text("Left", .5, .95)
    grid::grid.text("A", .9, .5) 
    grid::grid.text("P", .1, .5)
    grid::grid.text("S", .5, .9)
    grid::grid.text("I", .5, .1)
    dev.off()
  }
  if (any(views == "right")) {
    rgl::par3d(userMatrix = x90n %*% z90p) # right view
    rgl::rgl.snapshot(paste(statdir, "right.png", sep = ""), fmt = "png", top = TRUE)
    aa <- png::readPNG(paste(statdir, "right.png", sep = ""))
    png(paste(statdir, "right.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
    grid::grid.raster(aa)
    grid::grid.text("Right", .5, .95)
    grid::grid.text("P", .9, .5)
    grid::grid.text("A", .1, .5)
    grid::grid.text("S", .5, .9)
    grid::grid.text("I", .5, .1)
    dev.off()
  }
  if (any(views == "inferior")) {
    rgl::par3d(userMatrix = x180p) # inferior view
    rgl::rgl.snapshot(paste(statdir, "inferior.png", sep = ""), fmt = "png", top = TRUE)
    aa <- png::readPNG(paste(statdir, "inferior.png", sep = ""))
    png(paste(statdir, "inferior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
    grid::grid.raster(aa)
    grid::grid.text("Inferior", .5, .95)
    grid::grid.text("L", .9, .5) 
    grid::grid.text("R", .1, .5)
    grid::grid.text("A", .5, .9)
    grid::grid.text("P", .5, .1)
    dev.off()
  }
  if (any(views == "posterior")) {
    rgl::par3d(userMatrix = x90p %*% y180p) # posterior view
    rgl::rgl.snapshot(paste(statdir, "posterior.png", sep = ""), fmt = "png", top = TRUE)
    aa <- png::readPNG(paste(statdir, "posterior.png", sep = ""))
    png(paste(statdir, "posterior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
    grid::grid.raster(aa)
    grid::grid.text("Posterior", .5, .95)
    grid::grid.text("R", .9, .5) 
    grid::grid.text("L", .1, .5)
    grid::grid.text("S", .5, .9)
    grid::grid.text("I", .5, .1)
    dev.off()
  }
  if (any(views == "superior")) {
    rgl::par3d(userMatrix = z180p) # superior view
    rgl::rgl.snapshot(paste(statdir, "superior.png", sep = ""), fmt = "png", top = TRUE)
    aa <- png::readPNG(paste(statdir, "superior.png", sep = ""))
    png(paste(statdir, "superior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
    grid::grid.raster(aa)
    grid::grid.text("Superior", .5, .95)
    grid::grid.text("R", .9, .5) 
    grid::grid.text("L", .1, .5)
    grid::grid.text("A", .5, .9)
    grid::grid.text("P", .5, .1)
    dev.off()
  }
  # Create anterior view image
  
  #lim$usr[1] = left x
  #lim$usr[2] = y bottom
  #lim$usr[3] = right x
  #lim$usr[4] = y top
  
  if (nviews == 3 | nviews
  
  gridView1 <- png::readPNG(paste(statdir, views[1], ".png", sep = ""))
  for (i in 2:nviews) {
    view_name <- paste(statdir, views[i], ".png", sep = "")
    mypic <- png::readPNG(view_name)
    if (i ==2 ) {
      gridView1 <- abind::abind(gridView1, mypic, along = 1)
    } else if (i == 3) {
      gridView1 <- abind::abind(gridView1, mypic, along = 1)
    } else if (i == 4) {
      if (nviews == 4 | nviews == 8) {
        gridView1 <- abind::abind(gridView1, mypic, along = 1)
      } else {
        
      }
    } else if (
  
  
  abcdef <- abind::abind(abind::abind(abind::abind(ff, dd, along = 1),
                                      abind::abind(aa, ee, along = 1), along = 2),
                         abind::abind(bb, cc, along = 1), along = 2)
  
  png(paste(statdir, "FullView.png", sep = ""), width = dim(abcdef)[2], height = dim(abcdef)[1])
  grid::grid.raster(abcdef)
  dev.off()
}
