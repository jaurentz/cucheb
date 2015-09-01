#include <cuchebpoly.h>

/* routine for creating step filter */
int cuchebpoly_stepfilter(double a, double b, double c, double d, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // check c and d
  if ( c >= d ) {
    return 1;
  }

  // compute lower bound 
  double lb;
  if (c <= a) {lb = a;}
  else if (c >= b) {return 1;}
  else {lb = c;}

  // compute upper bound 
  double ub;
  if (d >= b) {ub = b;}
  else if (d <= a) {return 1;}
  else {ub = d;}

  // set degree
  ccp->degree = MAX_DOUBLE_DEG;

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // scale everything to [-1,1]
  lb = (2.0*lb - (b+a))/(b-a);
  ub = (2.0*ub - (b+a))/(b-a);
 
  // compute Chebyshev coefficients
  double pi = DOUBLE_PI;
  double aclb = acos(lb);
  double acub = acos(ub);
  ccp->coeffs[0] = (aclb - acub)/pi;
  for (int ii=1; ii<MAX_DOUBLE_DEG+1; ii++) {
    ccp->coeffs[ii] = 2.0*(sin(ii*aclb) - sin(ii*acub))/(ii*pi);
  }

  // apply Jackson damping
  int deg = MAX_DOUBLE_DEG;
  double alpha = 1.0/(deg+2.0);
  double beta = sin(pi*alpha);
  double gamma = cos(pi*alpha);
  for (int ii=0; ii<MAX_DOUBLE_DEG+1; ii++) {
    ccp->coeffs[ii] = alpha*((deg+2.0-ii)*beta*cos(ii*pi*alpha) +
                       sin(ii*pi*alpha)*gamma)*ccp->coeffs[ii]/beta;
  }

  // chop Chebyshev coefficients
  cuchebpoly_chop(&(ccp->degree),&(ccp->coeffs[0]));

  // return 
  return 0;

}

