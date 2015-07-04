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

  // compute function values for f(x) = 2*x^2 - 1
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    ccp1.coeffs[ii] = 2.0*pow(ccp1.coeffs[ii],2) - 1.0;
  }
  cuchebpoly_printlong(&ccp1);
 
  // return 
  return 0;

}
