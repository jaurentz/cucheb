#include <cucheb.h>

/* routine for Chebyshev points */
int cuchebpoly_points(double a, double b, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set points
  double alpha = (b-a)/2.0;
  double beta = (b+a)/2.0;
  double* points;
  points = &(ccp->points)[0];
  for (int ii=0; ii < 2*MAX_DOUBLE_DEG; ii++) {
    points[ii] = alpha*sin(DOUBLE_PI*(MAX_DOUBLE_DEG-2.0*ii)/(2.0*MAX_DOUBLE_DEG)) +
                  beta;
  }

  // return 
  return 0;

}

