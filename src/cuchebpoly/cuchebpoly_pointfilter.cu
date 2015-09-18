#include <cucheb.h>

/* routine for creating point filter */
int cuchebpoly_pointfilter(double a, double b, double rho, int order, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // compute shift
  double shift;
  if (rho <= a) {shift = a;}
  else if (rho >= b) {shift = b;}
  else {shift = rho;}

  // map shift to [-1,1]
  shift = (2.0*shift - b - a)/(b-a);

  // use clenshaw algorithm to compute coefficients for delta function
  double* coeffs;
  coeffs = &(ccp->coeffs)[0];
  coeffs[0] = 1.0;
  coeffs[1] = shift;
  for (int ii=2; ii < DOUBLE_DEG; ii++) {
    coeffs[ii] = 2.0*shift*coeffs[ii-1] - coeffs[ii-2];
  }
  coeffs[DOUBLE_DEG] = shift*coeffs[DOUBLE_DEG-1] - coeffs[DOUBLE_DEG-2];
 
  // set length
  ccp->degree = min(MAX_DOUBLE_DEG,max(0,order));

  // return 
  return 0;

}

