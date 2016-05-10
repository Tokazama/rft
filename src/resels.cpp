#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector resels(IntegerVector m, NumericVector fwhm, IntegerVector DIM) {
  int D = DIM.size();
  NumericVector resels(4);
  int v = m.size();
  
  int xm, ym, zm, xym, yzm, xz, xyzm;
  int Ex, Ey, Ez, Fxy, Fyz, Fxz, Fxyz;
  double rx = 1 / fwhm(0);
  double ry = 1 / fwhm(1);
  if (D == 2) {
    double rz = 0;
    for (int i = 0; i < v; i++ ) {
      if (m(i) == 1) {
        IntegerVector point = indx_helper(i, DIM(0), DIM(1));
        if ( (point(1) > 0) && (point(1) < v) && (m(point(1)) == 1) )
          Ex++;
        if ( (point(3) > 0) && (point(3) < v) && (m(point(3)) == 1) )
          Ey++;
        if ( ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(3) > 0) && (point(3) < v) && (m(point(3)) == 1)) &&
             ((point(5) > 0) && (point(5) < v) && (m(point(5)) == 1)) )
          Fxy++;
      }
    }
  } else if (D == 3) {
    double rz = 1 / fwhm(2);
    for (int i = 0; i < v; i++ ) {
      if (m(i) == 1) {
        IntegerVector point = indx_helper(i, DIM(0), DIM(1), DIM(2));
        if ( (point(1) > 0) && (point(1) < v) && (m(point(1)) == 1) )
          Ex++;
        if ( (point(3) > 0) && (point(3) < v) && (m(point(3)) == 1) )
          Ey++;
        if ( (point(5) > 0) && (point(5) < v) && (m(point(5)) == 1) )
          Ez++;
        if ( ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(3) > 0) && (point(3) < v) && (m(point(3)) == 1)) &&
             ((point(7) > 0) && (point(7) < v) && (m(point(7)) == 1)) )
          Fxy++;
        if ( ((point(3) > 0) && (point(3) < v) && (m(point(3)) == 1)) &&
             ((point(5) > 0) && (point(5) < v) && (m(point(5)) == 1)) &&
             ((point(11) > 0) && (point(11) < v) && (m(point(11)) == 1)) )
          Fyz++;
        if ( ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(5) > 0) && (point(5) < v) && (m(point(5)) == 1)) &&
             ((point(15) > 0) && (point(15) < v) && (m(point(15)) == 1)) )
          Fxz++;
        if ( ((point(1) > 0) && (point(1) < v) && (m(point(1)) == 1)) &&
             ((point(3) > 0) && (point(3) < v) && (m(point(3)) == 1)) &&
             ((point(5) > 0) && (point(5) < v) && (m(point(5)) == 1)) &&
             ((point(7) > 0) && (point(7) < v) && (m(point(7)) == 1)) &&
             ((point(11) > 0) && (point(11) < v) && (m(point(11)) == 1)) &&
             ((point(15) > 0) && (point(15) < v) && (m(point(15)) == 1)) &&
             ((point(19) > 0) && (point(19) < v) && (m(point(19)) == 1)) )
          Fxyz++;
      }
    }
  }
  resels(0) = v - (Ex +Ey + Ez) + (Fyz + Fxz + Fxy) - Fxyz;
  resels(1) = (Ex - Fxy - Fxz + Fxyz) * rx + (Ey - Fxy - Fyz + Fxyz) * ry + (Ez - Fxz - Fyz + Fxyz) * rz;
  resels(2) = (Fxy - Fxyz) * rx * ry + (Fxz - Fxyz) * rx * rz + (Fyz - Fxyz) * ry * rz;
  resels(3) = Fxyz * rx * ry * rz;
  return resels;
  }
