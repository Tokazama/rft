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
rg::par3d(userMatrix = x180p) # inferior view
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

png(paste(statdir, "anterior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
grid::grid.raster(aa)
title(main = "Anterior")
#lim <-par()
#xmidline <- mean(lim$usr[1:2])
#ymidline <- mean(lim$usr[3:4])
#text(xmidline, (lim$usr[3] + .01 * ymidline), "I") # inferior label
#plot(1:dim(aa)[1], 1:dim(aa)[2], type = "n", xaxt = 'n', yaxt = 'n', ann = FALSE, asp = 1)

bb <- png::readPNG(paste(statdir, "left.png", sep = ""))
png(paste(statdir, "left.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
grid::grid.raster(bb)
title(main = "Left")

cc <- png::readPNG(paste(statdir, "right.png", sep = ""))
png(paste(statdir, "right.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
grid::grid.raster(c)
title(main = "Right")

dd <- png::readPNG(paste(statdir, "inferior.png", sep = ""))
png(paste(statdir, "inferior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
grid::grid.raster(dd)
title(main = "Inferior")

ee <- png::readPNG(paste(statdir, "posterior.png", sep = ""))
png(paste(statdir, "posterior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
grid::grid.raster(ee)
title(main = "Posterior")

ff <- png::readPNG(paste(statdir, "superior.png", sep = ""))
png(paste(statdir, "superior.png", sep = ""), width = dim(aa)[2], height = dim(aa)[1])
grid::grid.raster(bb)
title(main = "Superior")

abcdef <- abind::abind(abind::abind(abind::abind(ff, dd, along = 1), # superior/inferior
                       abind::abind(aa, ee, along = 1), along = 2), # anterior/posterior
                       abind::abind(bb, cc, along = 1), along = 2) # left/right

png(paste(statdir, "FullView.png", sep = ""), width = dim(abcdef)[2], height = dim(abcdef)[1])
grid::grid.raster(abcdef)
dev.off()
}            

