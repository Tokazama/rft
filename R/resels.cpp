#include <RcppArmadillo.h>
using namespace arma ;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec resels(arma::cube mask, arma::vec fwhm) {
  
  int x2 = mask.n_rows ;
  int y2 = mask.n_cols ;
  int z2 = mask.n_slices ;
  int x1 = x2 - 1 ;
  int y1 = y2 - 1 ;
  int z1 = z2 - 1 ;
  
  double vox = accu(mask) ;
  double rx = mask.n_rows / fwhm(1) ;
  double ry = mask.n_cols / fwhm(2) ;
  double rz = mask.n_slices / fwhm(3) ;

  arma::cube xm = mask( span(0, x1), span(0, y1), span(0, z1)) +
                  mask( span(1, x2), span(0, y1), span(0, z1)) ;

  arma::cube ym = mask( span(0, x1), span(0, y1), span(0, z1)) +
                  mask( span(0, x1), span(1, y2), span(0, z1)) ;

  arma::cube zm = mask( span(0, x1), span(0, y1), span(0, z1)) +
                  mask( span(0, x1), span(0, y1), span(1, z2)) ;

  arma::cube xym = mask( span(0, x1), span(0, y1), span(0, z1)) +
                   mask( span(1, x2), span(0, y1), span(0, z1)) +
                   mask( span(0, x1), span(1, y2), span(0, z1)) +
                   mask( span(1, x2), span(1, y2), span(0, z1)) ;

  arma::cube xzm = mask( span(0, x1), span(0, y1), span(0, z1)) +
                   mask( span(1, x2), span(0, y1), span(0, z1)) +
                   mask( span(0, x1), span(0, y1), span(1, z2)) +
                   mask( span(1, x2), span(0, y1), span(1, z2)) ;

  arma::cube yzm = mask( span(0, x1), span(0, y1), span(0, z1)) +
                   mask( span(0, x1), span(1, y2), span(0, z1)) +
                   mask( span(0, x1), span(0, y1), span(1, z2)) +
                   mask( span(0, x1), span(1, y2), span(1, z2)) ;

  arma::cube xyzm = mask( span(0, x1), span(0, y1), span(0, z1)) +
                    mask( span(1, x2), span(0, y1), span(0, z1)) +
                    mask( span(0, x1), span(1, y2), span(0, z1)) +
                    mask( span(0, x1), span(0, y1), span(1, z2)) +
                    mask( span(1, x2), span(1, y2), span(0, z1)) +
                    mask( span(0, x1), span(1, y2), span(1, z2)) +
                    mask( span(1, x2), span(0, y1), span(1, z2)) +
                    mask( span(1, x2), span(1, y2), span(1, z2)) ;

  double Ex = accu((xm == 2)) / 2 ;
  double Ey = accu((ym == 2)) / 2 ;
  double Ez = accu((zm == 2)) / 2 ;
  double Fxy = accu((xym == 4)) / 4 ;
  double Fxz = accu((xzm == 4)) / 4 ;
  double Fyz = accu((yzm == 4)) / 4 ;
  double Fxyz = accu((xyzm == 8)) / 8 ;

  arma::vec resels(4) ;

  resels(1) = (vox - (Ex + Ey + Ez) + (Fyz + Fxz + Fxy) - Fxyz) ;
  resels(2) = (Ex - Fxy - Fxz + Fxyz) * rx + (Ey - Fxy - Fyz + Fxyz) * ry + (Ez - Fxz - Fyz + Fxyz)*rz ;
  resels(3) = (Fxy - Fxyz) * rx * ry + (Fxz - Fxyz) * rx * rz + (Fyz - Fxyz) * ry * rz ;
  resels(4) = Fxyz * rx * ry * rz ;
  
  return resels ;
}
