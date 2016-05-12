#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector resels(IntegerVector m, NumericVector fwhm, NumericVector DIM) {
  double D = DIM.size();
  NumericVector resels(4);
  double v = m.size();
  double Ex = 0.0, Ey = 0.0, Ez = 0.0, Fxy = 0.0, Fyz = 0.0, Fxz = 0.0, Fxyz = 0.0;
  double rx = 1.0 / fwhm(0);
  double ry = 1.0 / fwhm(1);
  double rz;
  NumericVector point;
  
  if (D == 2) {
    NumericVector point_indx(3);
    point_indx(0)  = (0.0 + 1.0) + ((0.0 + 0.0) * DIM(0));  // 322
    point_indx(1)  = (0.0 + 0.0) + ((0.0 + 1.0) * DIM(0));  // 232
    point_indx(2)  = (0.0 + 1.0) + ((0.0 + 1.0) * DIM(0));  // 332
    
    rz = 0.0;
    for (double i = 0.0; i < v; i++ ) {
      if (m(i) == 1) {
        point = i + point_indx;
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
    rz = 1.0 / fwhm(2);
    NumericVector point_indx(26);
    // x:
    point_indx(0)  = (0.0 - 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 122
    point_indx(1)  = (0.0 + 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 322
    // y:
    point_indx(2)  = (0.0 + 0.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 212
    point_indx(3)  = (0.0 + 0.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 232
    // z:
    point_indx(4)  = (0.0 + 0.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 221
    point_indx(5)  = (0.0 + 0.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 223
    // xy:
    point_indx(6)  = (0.0 - 1.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 112
    point_indx(7)  = (0.0 + 1.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 332
    point_indx(8)  = (0.0 - 1.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 132
    point_indx(9)  = (0.0 + 1.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 312
    // yz:
    point_indx(10) = (0.0 + 0.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 211
    point_indx(11) = (0.0 + 0.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 233
    point_indx(12) = (0.0 + 0.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 213
    point_indx(13) = (0.0 + 0.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 231
    // xz:
    point_indx(14) = (0.0 - 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 121
    point_indx(15) = (0.0 + 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 323
    point_indx(16) = (0.0 - 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 123
    point_indx(17) = (0.0 + 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 321
    // xyz:
    point_indx(18) = (0.0 - 1.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 111
    point_indx(19) = (0.0 + 1.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 333
    point_indx(20) = (0.0 - 1.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 133
    point_indx(21) = (0.0 + 1.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 311
    point_indx(22) = (0.0 - 1.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 113
    point_indx(23) = (0.0 + 1.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 331
    point_indx(24) = (0.0 - 1.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 - 1.0) * DIM(0) * DIM(1));  // 131
    point_indx(25) = (0.0 + 1.0) + ((0.0 - 1.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 313
    for (double i = 0; i < v; i++ ) {
      if (m(i) == 1) {
        point = i + point_indx;
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
