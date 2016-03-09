#' Estimates image resels
#' 
#' Utilize the estimated FWHM to find the resels per voxel
#' 
#' @param mask statistical value (typically the maxima of a cluster or SPM)
#' @param fwhm degrees of freedom expressed as df[degrees of interest, degrees of error]
#' @return A vector of resels for dimensions 0:D
#'
#' @seealso rftPval, euler
#' @note: function currently in beta phase
#' @example 
#' 
#' # generate some data as if we just fitted a linear regression
#' outimg1 <- makeImage(c(10, 10, 10), rt(1000))
#' mask <- getMask(outimg1)
#' fwhm <- estSmooth(outimg1, mask)
#' resels <- resels(maskimg, fwhm$fwhm)
#' ec <- euler(max(outimg1), c(1, 1) fieldType = "T")
#' pvox <- sum(ec * resels)
#'
#' @export resels
resels <- function(mask, fwhm, quick = TRUE) {
  dimx <- dim(mask)[1]
  dimy <- dim(mask)[2]
  dimz <- dim(mask)[3]
  mask <- iMath(mask, "PadImage", 1)
  if (quick == "TRUE") {
    .Call("resels", mask, fwhm)
  } else {
      nvox <- sum(as.array(mask))
      x1 <- 2:(dimx + 1)
      y1 <- 2:(dimy + 1)
      z1 <- 2:(dimz + 1)
      x2 <- 3:(dimx + 2)
      y2 <- 3:(dimy + 2)
      z2 <- 3:(dimz + 2)
      rx <- 1/ (fwhm[1])
      ry <- 1/ (fwhm[2])
      rz <- 1/ (fwhm[3])
      
      m <- mask[x1, y1, z1]
      xm   <- m + mask[x2, y1, z1]
      ym   <- m + mask[x1, y2, z1]
      zm   <- m + mask[x1, y1, z2]
      xym  <- m + mask[x2, y1, z1] + mask[x1, y2, z1] + mask[x2, y2, z1]
      xzm  <- m + mask[x2, y1, z1] + mask[x1, y1, z2] + mask[x2, y1, z2]
      yzm  <- m + mask[x1, y2, z1] + mask[x1, y1, z2] + mask[x1, y2, z2]
      xyzm <- m + mask[x2, y1, z1] + mask[x1, y2, z1] + mask[x1, y1, z2] + mask[x2, y2, z1] + mask[x2, y1, z2] + mask[x1, y2, z2] + mask[x2, y2, z2]
      
      # extract number of voxels that fits each set of parameters (see Worsley 1996 for exact definition of parameters)
      Ex   <- sum(xm[xm == 2]) / 2
      Ey   <- sum(ym[ym == 2]) / 2
      Ez   <- sum(zm[zm == 2]) / 2
      Fxy  <- sum(xym[xym == 4]) / 4
      Fxz  <- sum(xzm[xzm == 4]) / 4
      Fyz  <- sum(yzm[yzm == 4]) / 4
      Fxyz <- sum(xyzm[xyzm == 8]) / 8
      
      resels <- c(0, 0, 0, 0)
      resels[1] <- (nvox - (Ex + Ey + Ez) + (Fyz + Fxz + Fxy) - Fxyz)
      resels[2] <- (Ex - Fxy - Fxz + Fxyz) * rx + (Ey - Fxy - Fyz + Fxyz) * ry + (Ez - Fxz - Fyz + Fxyz) * rz
      resels[3] <- (Fxy - Fxyz) * rx * ry + (Fxz - Fxyz) * rx * rz + (Fyz - Fxyz) * ry * rz
      resels[4] <- Fxyz * rx * ry * rz
      resels
    }
  }
