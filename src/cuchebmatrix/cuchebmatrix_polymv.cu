#include <cuchebmatrix.h>

/* routine for poly mv multiply on GPU */
int cuchebmatrix_polymv(cuchebmatrix* ccm, cuchebpoly* ccp, double* alpha, double* x, 
                        double* beta, double* y){

  // local variables
  int n, deg;
  double a, b;
  double* coeffs;
  double* dvec1;
  double* dvec2;
  n = ccm->m;
  deg = ccp->degree;
  a = ccp->a;
  b = ccp->b;
  coeffs = &(ccp->coeffs)[0];
  dvec1 = &(ccm->dtemp)[0];
  dvec2 = &(ccm->dtemp)[n];

 
  // return 
  return 0;

}

