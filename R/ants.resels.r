ants.resels <- function(mask, fwhm) {

  P<-sum(as.array(mask))
  dimx <- dim(mask)[1]
  dimy <- dim(mask)[2]
  dimz <- dim(mask)[3]
  rx <- 1/(fwhm[1])
  ry <- 1/(fwhm[2])
  rz <- 1/(fwhm[3])

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
              if(mask[i,j,k]==1){
                  Ex <- ifelse(mask[i+1,j,k]==1,Ex+1,Ex)
                  Ey <- ifelse(mask[i,j+1,k]==1,Ey+1,Ey)
                  Ez <- ifelse(mask[i,j,k+1]==1,Ez+1,Ez)
                  Fxy <- ifelse(mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1,Fxy+1,Fxy)
                  Fxz <- ifelse(mask[i+1,j,k]==1 && mask[,j,k+1]==1 && mask[i+1,j,k+1]==1,Fxz+1,Fxz)
                  Fyz <- ifelse(mask[i,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i,j+1,k+1]==1,Fyz+1,Fyz)
                  cubes <- ifelse(mask[i,j,k]==1 && mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i+1,j,k+1]==1 && mask[i,j+1,k+1]==1 && mask[i+1,j+1,k+1]==1,cubes+1,cubes)
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
