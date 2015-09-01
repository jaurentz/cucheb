#include <cuchebpoly.h>

/* driver */
int main(){

  // cuhebpoly 1
  cuchebpoly ccp1;

  // initialize to identity and print
  cuchebpoly_init(&ccp1);
  cuchebpoly_print(&ccp1);
  cuchebpoly_printlong(&ccp1);

  // compute Chebyshev points in [-1,1]
  cuchebpoly_points(-1.0,1.0,&ccp1);
  cuchebpoly_printlong(&ccp1);

  // compute function values for f(x) = exp(-100*(x-.7)^2)
  for (int ii=0; ii<2*DOUBLE_DEG; ii++) {
    ccp1.points[ii] = exp(-100.0*pow(ccp1.points[ii]-.7,2));
  }
  cuchebpoly_printlong(&ccp1);
 
  // compute Chebyshev coefficients
  cuchebpoly_coeffs(&ccp1);
  cuchebpoly_printlong(&ccp1);

  // chop Chebyshev coefficients
  cuchebpoly_chop(&ccp1);
  cuchebpoly_print(&ccp1);

  // cuhebpoly 2
  cuchebpoly ccp2;

  // point_filter
  cuchebpoly_pointfilter(0.0,1000.0,900.0,&ccp2);
  cuchebpoly_print(&ccp2);

  // cuhebpoly 3
  cuchebpoly ccp3;

  // step_filter
  cuchebpoly_stepfilter(0.0,100.0,10.0,20.0,&ccp3);
  cuchebpoly_print(&ccp3);

  // free memory
  cuchebpoly_destroy(&ccp1);
  cuchebpoly_destroy(&ccp2);
  cuchebpoly_destroy(&ccp3);

  // return 
  return 0;

}
