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
  NumericVector Vxx(xv), Vyy(xv), Vzz(xv),
                Vxy(xv), Vyz(xv), Vxz(xv),
                rpv_img(xv), fwhm(3);
  double rpv;
  for ( int i = 0; i < xv; i++ ) {
    point = indx(i) + point_indx;
    for ( int j = 0; j < n; j++ ) {
      double dx = 0.0, dy = 0.0, dz = 0.0, dxy = 0.0, dyz = 0.0, dxz = 0.0, dxyz = 0.0, denom_x = 0.0, denom_y = 0.0,
        denom_z = 0.0, denom_xy = 0.0, denom_yz = 0.0, denom_xz = 0.0, denom_xyz = 0.0, dxx, dyy, dzz;
      if ( (point(0) >= 0) && (point(0) < v) && (m(point(0)) == 1) ) {
        dx = dx + std::abs(sr(j, indx(i)) - sr(j, point(0))) * 6.0;
        denom_x = denom_x + 6.0;
      }
      if ( (point(1) >= 0) && (point(1) < v) && (m(point(1)) == 1) ) {
        dx = dx + std::abs(sr(j, indx(i)) - sr(j, point(1))) * 6.0;
        denom_x = denom_x + 6.0;
      }
      if ( (point(2) >= 0) && (point(2) < v) && (m(point(2)) == 1) ) {
        dy = dy + std::abs(sr(j, indx(i)) - sr(j, point(2))) * 6.0;
        denom_y = denom_y + 6.0;
      }
      if ( (point(3) >= 0) && (point(3) < v) && (m(point(3)) == 1) ) {
        dy = dy + std::abs(sr(j, indx(i)) - sr(j, point(3))) * 6.0;
        denom_y = denom_y + 6.0;
      }
      if ( (point(4) >= 0) && (point(4) < v) && (m(point(4)) == 1) ) {
        dz = dz + std::abs(sr(j, indx(i)) - sr(j, point(4))) * 6.0;
        denom_z = denom_x + 6.0;
      }
      if ( (point(5) >= 0) && (point(5) < v) && (m(point(5)) == 1) ) {
        dz = dz + std::abs(sr(j, indx(i)) - sr(j, point(5))) * 6.0;
        denom_z = denom_z + 6.0;
      }
      if ( (point(6) >= 0) && (point(6) < v) && (m(point(6)) == 1) ) {
        dxy = dxy + std::abs(sr(j, indx(i)) - sr(j, point(6))) * 3.0;
        denom_xy = denom_xy + 3.0;
      }
      if ( (point(7) >= 0) && (point(7) < v) && (m(point(7)) == 1) ) {
        dxy = dxy + std::abs(sr(j, indx(i)) - sr(j, point(7))) * 3.0;
        denom_xy = denom_xy + 3.0;
      }
      if ( (point(8) >= 0) && (point(8) < v) && (m(point(8)) == 1) ) {
        dxy = dxy + std::abs(sr(j, indx(i)) - sr(j, point(8))) * 3.0;
        denom_xy = denom_xy + 3.0;
      }
      if ( (point(9) >= 0) && (point(9) < v) && (m(point(9)) == 1) ) {
        dxy = dxy + std::abs(sr(j, indx(i)) - sr(j, point(9))) * 3.0;
        denom_xy = denom_xy + 3.0;
      }
      if ( (point(10) >= 0) && (point(10) < v) && (m(point(10)) == 1) ) {
        dyz = dyz + std::abs(sr(j, indx(i)) - sr(j, point(10))) * 3.0;
        denom_yz = denom_xy + 3.0;
      }
      if ( (point(11) >= 0) && (point(11) < v) && (m(point(11)) == 1) ) {
        dyz = dyz + std::abs(sr(j, indx(i)) - sr(j, point(11))) * 3.0;
        denom_yz = denom_xy + 3.0;
      }
      if ( (point(12) >= 0) && (point(12) < v) && (m(point(12)) == 1) ) {
        dyz = dyz + std::abs(sr(j, indx(i)) - sr(j, point(12))) * 3.0;
        denom_yz = denom_xy + 3.0;
      }
      if ( (point(13) >= 0) && (point(13) < v) && (m(point(13)) == 1) ) {
        dyz = dyz + std::abs(sr(j, indx(i)) - sr(j, point(13))) * 3.0;
        denom_yz = denom_xy + 3.0;
      }
      if ( (point(14) >= 0) && (point(14) < v) && (m(point(14)) == 1) ) {
        dxz = dxz + std::abs(sr(j, indx(i)) - sr(j, point(14))) * 3.0;
        denom_xz = denom_xz + 3.0;
      }
      if ( (point(15) >= 0) && (point(15) < v) && (m(point(15)) == 1) ) {
        dxz = dxz + std::abs(sr(j, indx(i)) - sr(j, point(15))) * 3.0;
        denom_xz = denom_xz + 3.0;
      }
      if ( (point(16) >= 0) && (point(16) < v) && (m(point(16)) == 1) ) {
        dxz = dxz + std::abs(sr(j, indx(i)) - sr(j, point(16))) * 3.0;
        denom_xz = denom_xz + 3.0;
      }
      if ( (point(17) >= 0) && (point(17) < v) && (m(point(17)) == 1) ) {
        dxz = dxz + std::abs(sr(j, indx(i)) - sr(j, point(17))) * 3.0;
        denom_xz = denom_xz + 3.0;
      }
      if ( (point(18) >= 0) && (point(18) < v) && (m(point(18)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(18))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(19) >= 0) && (point(19) < v) && (m(point(19)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(19))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(20) >= 0) && (point(20) < v) && (m(point(20)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(20))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(21) >= 0) && (point(21) < v) && (m(point(21)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(21))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(22) >= 0) && (point(22) < v) && (m(point(22)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(22))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(23) >= 0) && (point(23) < xv) && (m(point(23)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(23))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(24) >= 0) && (point(24) < v) && (m(point(24)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(24))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      if ( (point(25) >= 0) && (point(25) < v) && (m(point(25)) == 1) ) {
        dxyz = dxyz + std::abs(sr(j, indx(i)) - sr(j, point(25))) * 2.0;
        denom_xyz = denom_xyz + 2.0;
      }
      dxx = (dx + dxy + dxz + dxyz); 
      dyy = (dy + dxy + dyz + dxyz); 
      dzz = (dz + dyz + dxz + dxyz); 
      if (dxx > 0)
        dxx = dxx / (denom_x + denom_xy + denom_xz + denom_xyz);
      if (dyy > 0)
        dyy = dyy / (denom_y + denom_xy + denom_yz + denom_xyz);
      if (dzz > 0)
        dzz = dzz / (denom_z + denom_yz + denom_xz + denom_xyz);
      
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
                       _["rpvImage"] = rpv_img);
}
