#include <cuchebpoly.h>

/* routine for Chebyshev points */
int cuchebpoints(double a, double b, double* points){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set points
  double alpha = (b-a)/2.0;
  double beta = (b+a)/2.0;
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    points[ii] = alpha*sin(DOUBLE_PI*(2.0*ii-DOUBLE_DEG)/(2.0*DOUBLE_DEG)) +
                  beta;
  }

  // return 
  return 0;

}

