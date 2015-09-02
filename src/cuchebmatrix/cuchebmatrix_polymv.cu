#include <cuchebmatrix.h>

/* routine for poly mv multiply on GPU */
/* y = p(A)*x */
int cuchebmatrix_polymv(cuchebmatrix* ccm, cuchebpoly* ccp, double* x, double* y){

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

  // scalars
  double zero = 0.0, mone = -1.0;
  double A, B;
  A = 4.0/(b-a);
  B = -2.0*(b+a)/(b-a);

  // initialize y
  cublasDcopy(ccm->cublashandle, n, x, 1, y, 1);
  cublasDscal(ccm->cublashandle, n, &coeffs[deg], y, 1);
 
  // initialize dvec1
  cublasDcopy(ccm->cublashandle, n, x, 1, dvec1, 1);
  cublasDscal(ccm->cublashandle, n, &zero, dvec1, 1);
 
  // loop for clenshaw
  for(int ii=0; ii<deg; ii++){

    // copy dvec1 to dvec2
    cublasDcopy(ccm->cublashandle, n, dvec1, 1, dvec2, 1);

    // copy y to dvec1
    cublasDcopy(ccm->cublashandle, n, y, 1, dvec1, 1);

    // scale A and B if ii == deg-1
    if(ii == deg-1){
      A = A/2.0;
      B = B/2.0;
    }

    // apply matrix
    cuchebmatrix_mv(ccm, &A, dvec1, &B, y);

    // add x
    cublasDaxpy(ccm->cublashandle, n, &coeffs[deg-ii-1], x, 1, y, 1);

    // subtract dvec2
    cublasDaxpy(ccm->cublashandle, n, &mone, dvec2, 1, y, 1);

  }

  // return 
  return 0;

}

