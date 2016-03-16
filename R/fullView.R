#'
#' bm <-getMask( antsImageRead( getANTsRData("ch2") ) )
#' brain <-renderSurfaceFunction( surfimg =list( bm ) , alphasurf=0.1, smoothsval = 1.5 )
#' fullView(tempfile( fileext='.png'))
#'
#' @export fullView

fullView <- function(statdir){
  
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
  rgl::par3d(userMatrix = x90n) # anterior view
  rgl::rgl.snapshot(paste(statdir, "anterior.png", sep = ""), fmt = "png",
                    top = TRUE)
  rgl::par3d(userMatrix = x90n %*% z90n) # left view
  rgl::rgl.snapshot(paste(statdir, "left.png", sep = ""), fmt = "png",
                    top = TRUE)
  rgl::par3d(userMatrix = x90n %*% z90p) # right view
  rgl::rgl.snapshot(paste(statdir, "right.png", sep = ""), fmt = "png",
                    top = TRUE)
  rgl::par3d(userMatrix = x180p) # inferior view
  rgl::rgl.snapshot(paste(statdir, "inferior.png", sep = ""), fmt = "png",
                    top = TRUE)
  rgl::par3d(userMatrix = x90p %*% y180p) # posterior view
  rgl::rgl.snapshot(paste(statdir, "posterior.png", sep = ""), fmt = "png",
                    top = TRUE)
  rgl::par3d(userMatrix = z180p) # superior view
  rgl::rgl.snapshot(paste(statdir, "superior.png", sep = ""), fmt = "png",
                    top = TRUE)
  
  # Create anterior view image
  aa <- png::readPNG(paste(statdir, "anterior.png", sep = ""))
  
  #lim$usr[1] = left x
  #lim$usr[2] = y bottom
  #lim$usr[3] = right x
  #lim$usr[4] = y top
  
  png(paste(statdir, "anterior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
  grid::grid.raster(aa)
  grid::grid.text("Anterior", .5, .95)
  grid::grid.text("L", .9, .5) 
  grid::grid.text("R", .1, .5)
  grid::grid.text("S", .5, .9)
  grid::grid.text("I", .5, .1)
  dev.off()
  
  bb <- png::readPNG(paste(statdir, "left.png", sep = ""))
  png(paste(statdir, "left.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
  grid::grid.raster(bb)
  grid::grid.text("Left", .5, .95)
  grid::grid.text("A", .9, .5) 
  grid::grid.text("P", .1, .5)
  grid::grid.text("S", .5, .9)
  grid::grid.text("I", .5, .1)
  dev.off()
  
  cc <- png::readPNG(paste(statdir, "right.png", sep = ""))
  png(paste(statdir, "right.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
  grid::grid.raster(cc)
  grid::grid.text("Right", .5, .95)
  grid::grid.text("P", .9, .5)
  grid::grid.text("A", .1, .5)
  grid::grid.text("S", .5, .9)
  grid::grid.text("I", .5, .1)
  dev.off()
  
  dd <- png::readPNG(paste(statdir, "inferior.png", sep = ""))
  png(paste(statdir, "inferior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
  grid::grid.raster(dd)
  grid::grid.text("Inferior", .5, .95)
  grid::grid.text("L", .9, .5) 
  grid::grid.text("R", .1, .5)
  grid::grid.text("A", .5, .9)
  grid::grid.text("P", .5, .1)
  dev.off()
  
  ee <- png::readPNG(paste(statdir, "posterior.png", sep = ""))
  png(paste(statdir, "posterior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
  grid::grid.raster(ee)
  grid::grid.text("Posterior", .5, .95)
  grid::grid.text("R", .9, .5) 
  grid::grid.text("L", .1, .5)
  grid::grid.text("S", .5, .9)
  grid::grid.text("I", .5, .1)
  dev.off()
  
  ff <- png::readPNG(paste(statdir, "superior.png", sep = ""))
  png(paste(statdir, "superior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
  grid::grid.raster(bb)
  grid::grid.text("Superior", .5, .95)
  grid::grid.text("R", .9, .5) 
  grid::grid.text("L", .1, .5)
  grid::grid.text("A", .5, .9)
  grid::grid.text("P", .5, .1)
  dev.off()
  
  aa <- png::readPNG(paste(statdir, "anterior.png", sep = ""))
  bb <- png::readPNG(paste(statdir, "left.png", sep = ""))
  cc <- png::readPNG(paste(statdir, "right.png", sep = ""))
  dd <- png::readPNG(paste(statdir, "inferior.png", sep = ""))
  ee <- png::readPNG(paste(statdir, "posterior.png", sep = ""))
  ff <- png::readPNG(paste(statdir, "superior.png", sep = ""))
  
  abcdef <- abind::abind(abind::abind(abind::abind(ff, dd, along = 1),
                                      abind::abind(aa, ee, along = 1), along = 2),
                         abind::abind(bb, cc, along = 1), along = 2)
  
  png(paste(statdir, "FullView.png", sep = ""), width = dim(abcdef)[2], height = dim(abcdef)[1])
  grid::grid.raster(abcdef)
  dev.off()
}
