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
  
  NumericVector point_indx(3);
  point_indx(0) = (0.0 + 1.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 322
  point_indx(1) = (0.0 + 0.0) + ((0.0 + 1.0) * DIM(0)) + ((0.0 + 0.0) * DIM(0) * DIM(1));  // 232
  point_indx(2) = (0.0 + 0.0) + ((0.0 + 0.0) * DIM(0)) + ((0.0 + 1.0) * DIM(0) * DIM(1));  // 223
  
  IntegerVector point(3);
  NumericVector Vxx(xv), Vyy(xv), Vzz(xv),  Vxy(xv), Vyz(xv), Vxz(xv), rpv_img(xv), fwhm(3);
  double dx, dy, dz, rpv, vtotal;
  for ( int i = 0; i < xv; i++ ) {
    point = indx(i) + point_indx;
    if ( (point(0) < v) && (point(0) >= 0) && (m(point(0)) == 1) &&
         (point(1) < v) && (point(1) >= 0) && (m(point(1)) == 1) &&
         (point(2) < v) && (point(2) >= 0) && (m(point(2)) == 1) ) {
      vtotal = vtotal + 1;
      for ( int j = 0; j < n; j++ ) {
        if ( (point(0) < v) && (point(0) >= 0) && (m(point(0)) == 1) )
          dx = sr(j, indx(i)) - sr(j, point(0));
        if ( (point(1) < v) && (point(1) >= 0) && (m(point(1)) == 1) )
          dy = sr(j, indx(i)) - sr(j, point(1));
        if ( (point(2) < v) && (point(2) >= 0) && (m(point(2)) == 1) )
          dz = sr(j, indx(i)) - sr(j, point(2));
        
        Vxx(i) = Vxx(i) + (dx * dx);
        Vyy(i) = Vyy(i) + (dy * dy);
        Vzz(i) = Vzz(i) + (dz * dz);
        Vxy(i) = Vxy(i) + (dx * dy);
        Vyz(i) = Vyz(i) + (dy * dz);
        Vxz(i) = Vxz(i) + (dx * dz);
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
    }
  }
  
  for ( int i = 0; i < xv; i++ ) {
    // estimate fwhm:
    fwhm(0) = fwhm(0) + (sqrt(Vxx(i) / (4.0 * log(2.0))) / vtotal);
    fwhm(1) = fwhm(1) + (sqrt(Vyy(i) / (4.0 * log(2.0))) / vtotal);
    fwhm(2) = fwhm(2) + (sqrt(Vzz(i) / (4.0 * log(2.0))) / vtotal);
  }
  rpv = pow((rpv / vtotal), (1.0 / 3.0));
  double fwhm_prod = (fwhm(0) * fwhm(1) * fwhm(2));
  fwhm(0) = rpv * (fwhm(0) / pow(fwhm_prod, (1.0 / 3.0)));
  fwhm(1) = rpv * (fwhm(1) / pow(fwhm_prod, (1.0 / 3.0)));
  fwhm(2) = rpv * (fwhm(2) / pow(fwhm_prod, (1.0 / 3.0)));
  fwhm = 1 / fwhm;
  return List::create( _["fwhm"] = fwhm,
                       _["rpvImage"] = rpv_img,
                       _["scale"] = scale,
                       _["rpv"] = rpv);
}
