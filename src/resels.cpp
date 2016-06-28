#include <Rcpp.h>
using namespace Rcpp;
  
// [[Rcpp::export]]
NumericVector resels(IntegerVector m, NumericVector fwhm, IntegerVector DIM) {
  int D = DIM.size();
  NumericVector resels(4);
  int dimx = DIM(0);
  int dimy = DIM(1);
  int dimz;
  if (D == 3)
    dimz = DIM(2);
  double v = m.size();
  double Ex = 0.0, Ey = 0.0, Ez = 0.0, Fxy = 0.0, Fyz = 0.0, Fxz = 0.0, Fxyz = 0.0;
  double rx = 1.0 / fwhm(0);
  double ry = 1.0 / fwhm(1);
  double rz;
  if (D == 3)
    rz = 1.0 / fwhm(2);
  double vtotal = 0.0;
  
  if (D == 2) {
    rz = 0.0;
    for (int i = 0.0; i < v; i++ ) {
      if (m(i) == 1) {
        int xc = (i % dimx);
        int yc = floor(i / dimx);
        yc = yc % dimy;
        
        IntegerVector point(7);
        point(0) = (xc + 1) + ((yc + 0) * DIM(0));
        point(1) = (xc + 0) + ((yc + 1) * DIM(0));
        point(2) = (xc + 1) + ((yc + 1) * DIM(0));
        
        vtotal++;
        if ( (point(0) > 0) && (point(0) < v) && (m(point(0)) == 1) )
          Ex++;
        if ( (point(1) > 0) && (point(1) < v) && (m(point(1)) == 1) )
          Ey++;
        if ( ((point(0) > 0) && (point(0) < v) && (m(point(0)) == 1)) &&
             ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(2) > 0) && (point(2) < v) && (m(point(2)) == 1)) )
          Fxy++;
      }
    }
  } else {
    for (int i = 0; i < v; i++ ) {
      if (m(i) == 1) {

        
        int xc = (i % dimx);
        
        int yc = floor(i / dimx);
        yc = yc % dimy;
        
        int zc = floor(i / (dimy * dimx));
        zc = zc % dimz;
        
        IntegerVector point(7);
        point(0) = (xc + 1) + ((yc + 0) * DIM(0)) + ((zc + 0) * DIM(0) * DIM(1));  // 322
        point(1) = (xc + 0) + ((yc + 1) * DIM(0)) + ((zc + 0) * DIM(0) * DIM(1));  // 232
        point(2) = (xc + 0) + ((yc + 0) * DIM(0)) + ((zc + 1) * DIM(0) * DIM(1));  // 223
        point(3) = (xc + 1) + ((yc + 1) * DIM(0)) + ((zc + 0) * DIM(0) * DIM(1));  // 332
        point(4) = (xc + 0) + ((yc + 1) * DIM(0)) + ((zc + 1) * DIM(0) * DIM(1));  // 233
        point(5) = (xc + 1) + ((yc + 0) * DIM(0)) + ((zc + 1) * DIM(0) * DIM(1));  // 323
        point(6) = (xc + 1) + ((yc + 1) * DIM(0)) + ((zc + 1) * DIM(0) * DIM(1));  // 333
        vtotal++;
        if ( (point(0) > 0) && (point(0) < v) && (m(point(0)) == 1) )
          Ex++;
        if ( (point(1) > 0) && (point(1) < v) && (m(point(1)) == 1) )
          Ey++;
        if ( (point(2) > 0) && (point(2) < v) && (m(point(2)) == 1) )
          Ez++;
        if ( ((point(0) > 0) && (point(0) < v) && (m(point(0)) == 1)) &&
             ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(3) > 0) && (point(3) < v) && (m(point(3)) == 1)) )
          Fxy++;
        if ( ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(2) > 0) && (point(2) < v) && (m(point(2)) == 1)) &&
             ((point(4) > 0) && (point(4) < v) && (m(point(4)) == 1)) )
          Fyz++;
        if ( ((point(0) > 0) && (point(0) < v) && (m(point(0)) == 1)) &&
             ((point(2) > 0) && (point(2) < v) && (m(point(2)) == 1)) &&
             ((point(5) > 0) && (point(5) < v) && (m(point(5)) == 1)) )
          Fxz++;
        if ( ((point(0) > 0) && (point(0) < v) && (m(point(0)) == 1)) &&
             ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(2) > 0) && (point(2) < v) && (m(point(2)) == 1)) &&
             ((point(3) > 0) && (point(3) < v) && (m(point(3)) == 1)) &&
             ((point(4) > 0) && (point(4) < v) && (m(point(4)) == 1)) &&
             ((point(5) > 0) && (point(5) < v) && (m(point(5)) == 1)) &&
             ((point(6) > 0) && (point(6) < v) && (m(point(6)) == 1)) )
          Fxyz++;
      }
    }
  }
  resels(0) = vtotal - (Ex + Ey + Ez) + (Fyz + Fxz + Fxy) - Fxyz;
  resels(1) = (Ex - Fxy - Fxz + Fxyz) * rx + (Ey - Fxy - Fyz + Fxyz) * ry + (Ez - Fxz - Fyz + Fxyz) * rz;
  resels(2) = (Fxy - Fxyz) * rx * ry + (Fxz - Fxyz) * rx * rz + (Fyz - Fxyz) * ry * rz;
  resels(3) = Fxyz * rx * ry * rz;
  return resels;
  }
