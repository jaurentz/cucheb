#include <cuchebpoly.h>

/* routine for chopping Chebyshev coefficients */
int cuchebchop(int* degree, double* coeffs){

  // find maximum
  double maximum = 0.0;
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    maximum = max(maximum, abs(coeffs[ii]));
  }

  // maximum == 0
  if (maximum == 0) {
    *degree = 0;
  }
  
  // maximum != 0
  else {
    
    // compute degree
    for (int ii=0; ii<MAX_DOUBLE_DEG+1; ii++) {

      // set current degree
      *degree = MAX_DOUBLE_DEG-ii;

      // exit if trailing coefficient is too large
      if (abs(coeffs[MAX_DOUBLE_DEG-ii]) >= DOUBLE_EPS*maximum) {
        break;
      }

    }

  }

  // return 
  return 0;

}

