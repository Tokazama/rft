#include <Rcpp.h>
using namespace Rcpp ;

// [[Rcpp::export]]
List estSmooth_test(NumericMatrix x, NumericVector m, double rdf, double nfull, NumericVector DIM, bool scaleResid) {
  // double D = DIM.size();
  double v = DIM(0) * DIM(1) * DIM(2);
  double n = x.nrow();
  double xv = x.ncol();
  double scale = (nfull / rdf) * (1.0 / n);
  // calculate standardized residuals and index active voxels:
  NumericMatrix sr(n, v);
  NumericVector indx(xv);
  double cntr = 0;
  for ( double i = 0; i < v; i++ ) {
    if ( m(i) == 1 ) {
      indx(cntr) = i;
      if (scaleResid == TRUE) {
        double mrss = pow(sum(pow(x( _, cntr), 2)), .5);
        sr( _, i) = x( _, cntr) / mrss;
      } else {
        sr( _, i) = x( _, cntr);
      }
      cntr = (cntr + 1);
    }
  }
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
  
  IntegerVector point(3);
  NumericVector Vxx(xv), Vyy(xv), Vzz(xv), Vxy(xv), Vyz(xv), Vxz(xv), rpv_img(xv), fwhm(3);
  double rpv;
  for ( int i = 0; i < xv; i++ ) {
    point = indx(i) + point_indx;
    for ( int j = 0; j < n; j++ ) {
      double dx1 = 0.0, dy1 = 0.0, dz1 = 0.0, denom_x = 0.0, denom_y = 0.0, denom_z = 0.0, denom_xy = 0.0, denom_yz = 0.0, denom_xz = 0.0, denom_xyz = 0.0;
      double dx2 = 0.0, dy2 = 0.0, dz2 = 0.0;
      double dxx = 0.0, dyy = 0.0, dzz = 0.0, dx = 0.0, dy = 0.0, dz = 0.0;
      // x:
      if ( (point(0) >= 0) && (point(0) < v) && (m(point(0)) == 1) )
        dx1 = (sr(j, indx(i)) - sr(j, point(0)));
      if ( (point(1) >= 0) && (point(1) < v) && (m(point(1)) == 1) )
        dx2 = (sr(j, indx(i)) - sr(j, point(1)));
      if ( (dx1 != 0) && (dx2 != 0) ) {
        dx = ((dx1 - dx2) / 2) * 6;
        denom_x = denom_x + 6;
      } else if ( (dx1 != 0) | (dx2 != 0) ) {
        dx = (dx1 - dx2) * 6; // only one will not equal zero
        denom_x = denom_x + 6;
      }
      
      // y:
      if ( (point(2) >= 0) && (point(2) < v) && (m(point(2)) == 1) )
        dy1 = (sr(j, indx(i)) - sr(j, point(2)));
      if ( (point(3) >= 0) && (point(3) < v) && (m(point(3)) == 1) )
        dy2 = (sr(j, indx(i)) - sr(j, point(3)));
      if ( (dy1 != 0) && (dy2 != 0) ) {
        dy = ((dy1 - dy2) / 2) * 6;
        denom_y = denom_y + 6;
      } else if ( (dx1 != 0) | (dx2 != 0) ) {
        dyy = dxx + (dy1 - dy2) * 6;
        denom_y = denom_y + 6;
      }
      
      // z:
      if ( (point(4) >= 0) && (point(4) < v) && (m(point(4)) == 1) )
        dz1 = (sr(j, indx(i)) - sr(j, point(4)));
      if ( (point(5) >= 0) && (point(5) < v) && (m(point(5)) == 1) )
        dz2 = (sr(j, indx(i)) - sr(j, point(5)));
      if ( (dz1 != 0) && (dz2 != 0) ) {
        dz = ((dz1 - dz2) / 2) * 6;
        denom_z = denom_z + 6;
      } else if ( (dz1 != 0) | (dz2 != 0) ) {
        dz = (dz1 - dz2) * 6;
        denom_z = denom_z + 6;
      }
      
      // xy:
      double dxy11 = 0.0, dxy12 = 0.0, dxy21 = 0.0, dxy22 = 0.0, dxy = 0.0;
      if ( (point(6) >= 0) && (point(6) < v) && (m(point(6)) == 1) )
        dxy11 = (sr(j, indx(i)) - sr(j, point(6)));
      if ( (point(7) >= 0) && (point(7) < v) && (m(point(7)) == 1) )
        dxy12 = (sr(j, indx(i)) - sr(j, point(7)));
      if ( (dxy11 != 0) && (dxy12 != 0) ) {
        dxy = dxy + ((dxy11 - dxy21) / 2) * 3;
        denom_x = denom_x + 3;
      } else if ( (dxy11 != 0) | (dxy12 != 0) ) {
        dxy = dxy + (dxy11 - dxy12) * 3;
        denom_xy = denom_xy + 3;
      }
      
      if ( (point(8) >= 0) && (point(8) < v) && (m(point(8)) == 1) )
        dxy21 = (sr(j, indx(i)) - sr(j, point(8)));
      if ( (point(9) >= 0) && (point(9) < v) && (m(point(9)) == 1) )
        dxy22 = (sr(j, indx(i)) - sr(j, point(9)));
      if ( (dxy21 != 0) && (dxy22 != 0) ) {
        dxy = dxy + ((dxy21 - dxy21) / 2) * 3;
        denom_x = denom_x + 3;
      } else if ( (dxy21 != 0) | (dxy22 != 0) ) {
        dxy = dxy + (dxy21 - dxy22) * 3;
        denom_xy = denom_xy + 3;
      }
      
      // yz:
      double dyz11 = 0.0, dyz12 = 0.0, dyz21 = 0.0, dyz22 = 0.0, dyz = 0.0;
      if ( (point(10) >= 0) && (point(10) < v) && (m(point(10)) == 1) )
        dyz11 = (sr(j, indx(i)) - sr(j, point(10)));
      if ( (point(11) >= 0) && (point(11) < v) && (m(point(11)) == 1) )
        dyz12 = (sr(j, indx(i)) - sr(j, point(11)));
      if ( (dyz11 != 0) && (dyz12 != 0) ) {
        dyz = dyz + ((dyz11 - dyz12) / 2) * 3;
        denom_yz = denom_yz + 3;
      } else if ( (dyz11 != 0) | (dyz12 != 0) ) {
        dyz = dyz + (dyz11 - dyz12) * 3;
        denom_yz = denom_yz + 3;
      }
      
      if ( (point(12) >= 0) && (point(12) < v) && (m(point(12)) == 1) )
        dyz21 = (sr(j, indx(i)) - sr(j, point(12)));
      if ( (point(13) >= 0) && (point(13) < v) && (m(point(13)) == 1) )
        dyz22 = (sr(j, indx(i)) - sr(j, point(13)));
      if ( (dyz21 != 0) && (dyz22 != 0) ) {
        dyz = dyz + ((dyz21 - dyz22) / 2) * 3;
        denom_yz = denom_yz + 3;
      } else if ( (dyz21 != 0) | (dxy22 != 0) ) {
        dyz = dyz + (dyz21 - dyz22) * 3;
        denom_yz = denom_yz + 3;
      }
      
      // xz:
      double dxz11 = 0.0, dxz12 = 0.0, dxz21 = 0.0, dxz22 = 0.0, dxz = 0.0;
      if ( (point(14) >= 0) && (point(14) < v) && (m(point(14)) == 1) )
        dxz11 = (sr(j, indx(i)) - sr(j, point(14)));
      if ( (point(15) >= 0) && (point(15) < v) && (m(point(15)) == 1) )
        dxz12 = (sr(j, indx(i)) - sr(j, point(15)));
      if ( (dxz11 != 0) && (dxz12 != 0) ) {
        dxz = dxz + ((dxz11 - dxz12) / 2) * 3;
        denom_xz = denom_xz + 3;
      } else if ( (dxz11 != 0) | (dxz12 != 0) ) {
        dxz = dxz + (dxz11 - dxz12) * 3;
        denom_xz = denom_xz + 3;
      }
      
      if ( (point(16) >= 0) && (point(16) < v) && (m(point(16)) == 1) )
        dxz21 = (sr(j, indx(i)) - sr(j, point(16)));
      if ( (point(17) >= 0) && (point(17) < v) && (m(point(17)) == 1) )
        dxz22 = (sr(j, indx(i)) - sr(j, point(17)));
      if ( (dxz21 != 0) && (dxz22 != 0) ) {
        dxz = dxz + ((dxz21 - dxz22) / 2) * 3;
        denom_xz = denom_xz + 3;
      } else if ( (dxz21 != 0) | (dxz22 != 0) ) {
        dxz = dxz + (dxz21 - dxz12) * 3;
        denom_xz = denom_xz + 3;
      }
      
      // xyz:
      double dxyz = 0.0, dxyz11 = 0.0, dxyz12 = 0.0, dxyz21 = 0.0, dxyz22 = 0.0, dxyz31 = 0.0, dxyz32 = 0.0, dxyz41 = 0.0, dxyz42 = 0.0;
      if ( (point(18) >= 0) && (point(18) < v) && (m(point(18)) == 1) )
        dxyz11 = (sr(j, indx(i)) - sr(j, point(18)));
      if ( (point(19) >= 0) && (point(19) < v) && (m(point(19)) == 1) )
        dxyz12 = (sr(j, indx(i)) - sr(j, point(19)));
      if ( (dxyz11 != 0) && (dxyz12 != 0) ) {
        dxyz = dxyz + ((dxyz11 - dxyz12) / 2) * 2;
        denom_xyz = denom_xyz + 2;
      } else if ( (dxyz11 != 0) | (dxyz12 != 0) ) {
        dxyz = dxyz + (dxyz11 - dxyz12) * 2;
        denom_xyz = denom_xyz + 2;
      }

      
      if ( (point(20) >= 0) && (point(20) < v) && (m(point(20)) == 1) )
        dxyz21 = (sr(j, indx(i)) - sr(j, point(20)));
      if ( (point(21) >= 0) && (point(21) < v) && (m(point(21)) == 1) )
        dxyz22 = (sr(j, indx(i)) - sr(j, point(21)));
      if ( (dxyz21 != 0) && (dxyz22 != 0) ) {
        dxyz = dxyz + ((dxyz21 - dxyz22) / 2) * 2;
        denom_xyz = denom_xyz + 2;
      } else if ( (dxyz21 != 0) | (dxyz22 != 0) ) {
        dxyz = dxyz + (dxyz21 - dxyz22) * 2;
        denom_xyz = denom_xyz + 2;
      }
      
      if ( (point(22) >= 0) && (point(22) < v) && (m(point(22)) == 1) )
        dxyz31 = (sr(j, indx(i)) - sr(j, point(22)));
      if ( (point(23) >= 0) && (point(23) < v) && (m(point(23)) == 1) )
        dxyz32 = (sr(j, indx(i)) - sr(j, point(23)));
      if ( (dxyz31 != 0) && (dxyz32 != 0) ) {
        dxyz = dxyz + ((dxyz31 - dxyz32) / 2) * 2;
        denom_xyz = denom_xyz + 2;
      } else if ( (dxyz31 != 0) | (dxyz32 != 0) ) {
        dxyz = dxyz + (dxyz31 - dxyz32) * 2;
        denom_xyz = denom_xyz + 2;
      }
      
      if ( (point(24) >= 0) && (point(24) < v) && (m(point(24)) == 1) )
        dxyz41 = (sr(j, indx(i)) - sr(j, point(24)));
      if ( (point(25) >= 0) && (point(25) < v) && (m(point(25)) == 1) )
        dxyz42 = (sr(j, indx(i)) - sr(j, point(25)));
      if ( (dxyz41 != 0) && (dxyz42 != 0) ) {
        dxyz = dxyz + ((dxyz41 - dxyz42) / 2) * 2;
        denom_xyz = denom_xyz + 2;
      } else if ( (dxyz41 != 0) | (dxyz42 != 0) ) {
        dxyz = dxyz + (dxyz41 - dxyz42) * 2;
        denom_xyz = denom_xyz + 2;
      }
      
      dxx = ( dx + dxy + dxz + dxyz ) / ( denom_x + denom_xy + denom_xz + denom_xyz );
      dyy = ( dy + dyz + dxy + dxyz ) / ( denom_y + denom_xy + denom_yz + denom_xyz );
      dzz = ( dz + dxz + dyz + dxyz ) / ( denom_z + denom_yz + denom_xz + denom_xyz );
      
      Vxx(i) = Vxx(i) + (dxx * dxx);
      Vyy(i) = Vyy(i) + (dyy * dyy);
      Vzz(i) = Vzz(i) + (dzz * dzz);
      Vxy(i) = Vxy(i) + (dxx * dyy);
      Vyz(i) = Vyz(i) + (dyy * dzz);
      Vxz(i) = Vxz(i) + (dxx * dzz);
    }
      Vxx(i) = Vxx(i) * scale;
      Vyy(i) = Vyy(i) * scale;
      Vzz(i) = Vzz(i) * scale;
      Vxy(i) = Vxy(i) * scale;
      Vyz(i) = Vyz(i) * scale;
      Vxz(i) = Vxz(i) * scale;
      
      rpv_img(i) = (Vxx(i) * Vyy(i) * Vzz(i)) +
                   (Vxy(i) * Vyz(i) * Vxz(i) * 2.0) -
                   (Vyz(i) * Vyz(i) * Vxx(i)) -
                   (Vxy(i) * Vxy(i) * Vzz(i)) -
                   (Vxz(i) * Vxz(i) * Vyy(i));
      if (rpv_img(i) < 0.0)
        rpv_img(i) = 0.0;
      else {
        rpv_img(i) = sqrt(rpv_img(i) / pow(4.0 * log(2.0), 3.0));
        rpv = rpv + rpv_img(i);
      }
      fwhm(0) = fwhm(0) + sqrt(Vxx(i) / (4.0 * log(2.0)));
      fwhm(1) = fwhm(1) + sqrt(Vyy(i) / (4.0 * log(2.0)));
      fwhm(2) = fwhm(2) + sqrt(Vzz(i) / (4.0 * log(2.0)));
    }
  fwhm = fwhm / xv;
  rpv = pow((rpv / xv), (1.0 / 3.0));
  double fwhm_prod = (fwhm(0) * fwhm(1) * fwhm(2));
  fwhm(0) = rpv * (fwhm(0) / pow(fwhm_prod, (1.0 / 3.0)));
  fwhm(1) = rpv * (fwhm(1) / pow(fwhm_prod, (1.0 / 3.0)));
  fwhm(2) = rpv * (fwhm(2) / pow(fwhm_prod, (1.0 / 3.0)));
  fwhm = 1 / fwhm;
  return List::create( _["fwhm"] = fwhm,
                       _["rpvImage"] = rpv_img,
                       _["rpv"] = rpv,
                       _["Vxx"] = Vxx,
                       _["Vyy"] = Vyy,
                       _["Vzz"] = Vzz);
}
