#include <cuchebpoly.h>

/* routine for basic initialization */
int cuchebpoly_init(cuchebpoly* ccp){

  // set degree
  ccp->degree = 0;
  
  // set a and b
  ccp->a = -1.0;
  ccp->b = 1.0;
  
  // set coeffs
  ccp->coeffs[0] = 1.0;
  for (int ii=0; ii<DOUBLE_DEG; ii++) {
    ccp->coeffs[ii+1] = 0.0;
  }

  // return 
  return 0;

}

