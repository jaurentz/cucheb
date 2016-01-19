#include <cucheb.h>

/* routine for chopping Chebyshev coefficients */
int cuchebpoly_chop(cuchebpoly* ccp){

  // find maximum
  double maximum = 0.0;
  for (int ii=0; ii<MAX_DOUBLE_DEG+1; ii++) {
    maximum = max(maximum, abs((ccp->coeffs)[ii]));
  }

  // maximum == 0
  if (maximum == 0) {
    ccp->degree = 0;
  }
  
  // maximum != 0
  else {
    
    // compute degree
    for (int ii=0; ii<MAX_DOUBLE_DEG+1; ii++) {

      // set current degree
      ccp->degree = MAX_DOUBLE_DEG-ii;

      // exit if trailing coefficient is too large
      if (abs((ccp->coeffs)[MAX_DOUBLE_DEG-ii]) >= DOUBLE_EPS*maximum) {
        break;
      }

    }

  }

  // return 
  return 0;

}

