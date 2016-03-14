#'
#' bm <-getMask( antsImageRead( getANTsRData("ch2") ) )
#' brain <-renderSurfaceFunction( surfimg =list( bm ) , alphasurf=0.1 ,smoothsval = 1.5 )
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



rgl::par3d(userMatrix = x90n) # anterior view
rgl::rgl.snapshot(paste(statdir, "anterior.png", sep = ""), fmt = "png",
                  top = TRUE)
rgl::par3d(userMatrix = x90n %*% z90n) # left view
rgl::rgl.snapshot(paste(statdir, "left.png", sep = ""), fmt = "png",
                  top = TRUE)
rgl::par3d(userMatrix = x90n %*% z90p) # right view
rgl::rgl.snapshot(paste(statdir, "right.png", sep = ""), fmt = "png",
                  top = TRUE)
rg::par3d(userMatrix = x180p) # inferior view
rgl::rgl.snapshot(paste(statdir, "inferior.png", sep = ""), fmt = "png",
                  top = TRUE)
rgl::par3d(userMatrix = x90p %*% y180p) # posterior view
rgl::rgl.snapshot(paste(statdir, "posterior.png", sep = ""), fmt = "png",
                  top = TRUE)
rgl::par3d(userMatrix = z180p) # superior view
rgl::rgl.snapshot(paste(statdir, "superior.png", sep = ""), fmt = "png",
                  top = TRUE)

aa <- png::readPNG(paste(statdir, "anterior.png", sep = ""))
bb <- png::readPNG(paste(statdir, "left.png", sep = ""))
cc <- png::readPNG(paste(statdir, "right.png", sep = ""))
dd <- png::readPNG(paste(statdir, "inferior.png", sep = ""))
ee <- png::readPNG(paste(statdir, "posterior.png", sep = ""))
ff <- png::readPNG(paste(statdir, "superior.png", sep = ""))

abcdef <- abind::abind(abind::abind(abind::abind(ff, dd, along = 1), # superior/inferior
                       abind::abind(aa, ee, along = 1), along = 2), # anterior/posterior
                       abind::abind(bb, cc, along = 1), along = 2) # left/right

png(paste(statdir, "FullView.png", sep = ""), width = dim(abcdef)[2], height = dim(abcdef)[1])
grid::grid.raster(abcdef)
dev.off()
}            

