#include <cuchebpoly.h>

/* routine for creating point filter */
int cuchebpoly_pointfilter(double a, double b, double rho, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // compute Chebyshev points in [a,b]
  cuchebpoly_points(a,b,&(ccp->coeffs[0]));

  // compute shift
  double shift;
  if (rho <= a) {shift = a;}
  else if (rho >= b) {shift = b;}
  else {shift = rho;}

  // compute function values for f(x) = exp(-100*(x-shift)^2)
  double scl = pow(b - a,2);
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    ccp->coeffs[ii] = exp(-100.0*pow(ccp->coeffs[ii]-shift,2)/scl);
  }
 
  // compute Chebyshev coefficients
  cuchebpoly_coeffs(&(ccp->coeffs[0]));

  // chop Chebyshev coefficients
  cuchebpoly_chop(&(ccp->degree),&(ccp->coeffs[0]));

  // return 
  return 0;

}

