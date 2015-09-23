#include <cucheb.h>

/* routine for creating gaussian filter */
int cuchebpoly_gaussianfilter(double a, double b, double rho, double tau, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // compute Chebyshev points in [a,b]
  cuchebpoly_points(a,b,ccp);

  // compute shift
  double shift;
  if (rho <= a) {shift = a;}
  else if (rho >= b) {shift = b;}
  else {shift = rho;}

  // compute function values for f(x) = exp(-100*(x-shift)^2)
  double scl = pow(b - a,2);
  for (int ii=0; ii < 2*MAX_DOUBLE_DEG; ii++) {
    (ccp->points)[ii] = exp(-abs(tau)*pow((ccp->points)[ii]-shift,2)/scl);
  }
 
  // compute Chebyshev coefficients
  cuchebpoly_coeffs(ccp);

  // chop Chebyshev coefficients
  cuchebpoly_chop(ccp);

  // return 
  return 0;

}

