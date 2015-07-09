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
  cuchebpoints(-1.0,1.0,&ccp1.coeffs[0]);
  cuchebpoly_printlong(&ccp1);

  // compute function values for f(x) = exp(-100*(x-.7)^2)
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    ccp1.coeffs[ii] = exp(-100.0*pow(ccp1.coeffs[ii]-.7,2));
  }
  cuchebpoly_printlong(&ccp1);
 
  // compute Chebyshev coefficients
  cuchebcoeffs(&ccp1.coeffs[0]);
  cuchebpoly_printlong(&ccp1);

  // chop Chebyshev coefficients
  cuchebchop(&ccp1.degree,&ccp1.coeffs[0]);
  cuchebpoly_print(&ccp1);

  // return 
  return 0;

}
