ants.resels <- function(mask, fwhm) {

  P<-labelGeometryMeasures(mask)[2,2]
  dimx <- dim(mask)[1]
  dimy <- dim(mask)[2]
  dimz <- dim(mask)[3]
  rx <- (dimx)/(fwhm[1])
  ry <- (dimy)/(fwhm[2])
  rz <- (dimz)/(fwhm[3])

  Ex <- 0
  Ey <- 0
  Ez <- 0
  Fxy <- 0
  Fxz <- 0
  Fyz <- 0
  cubes <- 0

  for (i in 1:(dimx)){
      for (j in 1:(dimy)){
          for (k in 1:(dimz)){
              if(getPixels(mask,i,j,k)[1]==1){
                  Ex <- ifelse(getPixels(mask,i+1,j,k)[1]==1,Ex+1,Ex)
                  Ey <- ifelse(getPixels(mask,i,j+1,k)[1]==1,Ey+1,Ey)
                  Ez <- ifelse(getPixels(mask,i,j,k+1)[1]==1,Ez+1,Ez)
                  Fxy <- ifelse(getPixels(mask,i+1,j,k)[1]==1 && getPixels(mask,i,j+1,k)[1]==1 && getPixels(mask,i+1,j+1,k)[1]==1,Fxy+1,Fxy)
                  Fxz <- ifelse(getPixels(mask,i+1,j,k)[1]==1 && getPixels(mask,i,j,k+1)[1]==1 && getPixels(mask,i+1,j,k+1)[1]==1,Fxz+1,Fxz)
                  Fyz <- ifelse(getPixels(mask,i,j+1,k)[1]==1 && getPixels(mask,i,j,k+1)[1]==1 && getPixels(mask,i,j+1,k+1)[1]==1,Fyz+1,Fyz)
                  cubes <- ifelse(getPixels(mask,i,j,k)[1]==1 && getPixels(mask,i+1,j,k)[1]==1 && getPixels(mask,i,j+1,k)[1]==1 && getPixels(mask,i+1,j+1,k)[1]==1 && getPixels(mask,i,j,k+1)[1]==1 && getPixels(mask,i+1,j,k+1)[1]==1 && getPixels(mask,i,j+1,k+1)[1]==1 && getPixels(mask,i+1,j+1,k+1)[1]==1,cubes+1,cubes)
              }
          }
      }
  }

  r1 <- (P-(Ex+Ey+Ez)+(Fyz+Fxz+Fxy)-cubes)
  r2 <- (((Ex-Fxy-Fxz+cubes)*rx)+((Ey-Fxy-Fyz+cubes)*ry)+((Ez-Fxz-Fyz+cubes)*rz))
  r3 <- (((Fxy-cubes)*rx*ry)+((Fxz-cubes)*rx*rz)+((Fyz-cubes)*ry*rz))
  r4 <- (cubes*rx*ry*rz)
  resel<-c(r1,r2,r3,r4)
  return(resel)
}
