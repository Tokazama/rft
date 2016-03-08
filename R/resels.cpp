#include <RcppArmadillo.h>
using namespace arma ;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::vec resels(arma::cube mask, arma::vec fwhm){
  int dimx = mask.n_rows ;
  int dimy = mask.n_cols ;
  int dimz = mask.n_slices ;
  double vox = accu(mask) ;
  double rx = mask.n_rows / fwhm(1) ;
  double ry = mask.n_cols / fwhm(2) ;
  double rz = mask.n_slices / fwhm(3) ;
  
  arma::cube xm = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                  mask( span(3, dimx + 2), span(2, dimy + 1), span(2, dimz + 1)) ;
  
  arma::cube ym = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                  mask( span(2, dimx + 1), span(3, dimy + 2), span(2, dimz + 1)) ;
  
  arma::cube zm = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                  mask( span(2, dimx + 1), span(2, dimy + 1), span(3, dimz + 2)) ;
  
  arma::cube xym = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                   mask( span(3, dimx + 2), span(2, dimy + 1), span(2, dimz + 1)) +
                   mask( span(2, dimx + 1), span(3, dimy + 2), span(2, dimz + 1)) +
                   mask( span(3, dimx + 2), span(3, dimy + 2), span(2, dimz + 1)) ;
                
  arma::cube xzm = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                   mask( span(3, dimx + 2), span(2, dimy + 1), span(2, dimz + 1)) +
                   mask( span(2, dimx + 1), span(2, dimy + 1), span(3, dimz + 2)) +
                   mask( span(3, dimx + 2), span(2, dimy + 1), span(3, dimz + 2)) ;
  
  arma::cube yzm = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                   mask( span(2, dimx + 1), span(3, dimy + 2), span(2, dimz + 1)) +
                   mask( span(2, dimx + 1), span(2, dimy + 1), span(3, dimz + 2)) +
                   mask( span(2, dimx + 1), span(3, dimy + 2), span(3, dimz + 2)) ;
  
  arma::cube xyzm = mask( span(2, dimx + 1), span(2, dimy + 1), span(2, dimz + 1)) +
                    mask( span(3, dimx + 2), span(2, dimy + 1), span(2, dimz + 1)) +
                    mask( span(2, dimx + 1), span(3, dimy + 2), span(2, dimz + 1)) +
                    mask( span(2, dimx + 1), span(2, dimy + 1), span(3, dimz + 2)) +
                    mask( span(3, dimx + 2), span(3, dimy + 2), span(2, dimz + 1)) +
                    mask( span(2, dimx + 1), span(3, dimy + 2), span(3, dimz + 2)) +
                    mask( span(3, dimx + 2), span(2, dimy + 1), span(3, dimz + 2)) +
                    mask( span(3, dimx + 2), span(3, dimy + 2), span(3, dimz + 2)) ;
  
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
