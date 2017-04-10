#include <cucheb.h>

/* routine for poly mm multiply on GPU */
/* Y = p(A)*X */
int cuchebmatrix_polymm(cuchebmatrix* ccm, cuchebpoly* ccp, int bsize,
                        double* X, double* Y, double* V1, double* V2){

  // local variables
  int n, deg;
  double a, b;
  double* coeffs;
  n = ccm->m;
  deg = ccp->degree;
  a = ccp->a;
  b = ccp->b;
  coeffs = &(ccp->coeffs)[0];

  // scalars
  double zero = 0.0, mone = -1.0;
  double A, B;
  A = 4.0/(b-a);
  B = -2.0*(b+a)/(b-a);

  // initialize Y
  cublasDcopy(ccm->cublashandle, bsize*n, X, 1, Y, 1);
  cublasDscal(ccm->cublashandle, bsize*n, &coeffs[deg], Y, 1);
 
  // initialize V1
  cublasDcopy(ccm->cublashandle, bsize*n, X, 1, V1, 1);
  cublasDscal(ccm->cublashandle, bsize*n, &zero, V1, 1);
 
clock_t tick;
double cp_time = 0.0;
double mm_time = 0.0;
double ad_time = 0.0;

  // loop for clenshaw
  for(int ii=0; ii<deg; ii++){

tick = clock();
    // copy V1 to V2
    cublasDcopy(ccm->cublashandle, bsize*n, V1, 1, V2, 1);

    // copy Y to V1
    cublasDcopy(ccm->cublashandle, bsize*n, Y, 1, V1, 1);
cp_time += (clock()-tick)/((double)CLOCKS_PER_SEC);

tick = clock();
    // scale A and B if ii == deg-1
    if(ii == deg-1){
      A = A/2.0;
      B = B/2.0;
    }

    // apply matrix
    cuchebmatrix_mm(ccm, bsize, &A, V1, &B, Y);
mm_time += (clock()-tick)/((double)CLOCKS_PER_SEC);

tick = clock();
    // add x
    cublasDaxpy(ccm->cublashandle, bsize*n, &coeffs[deg-ii-1], X, 1, Y, 1);

    // subtract V2
    cublasDaxpy(ccm->cublashandle, bsize*n, &mone, V2, 1, Y, 1);
ad_time += (clock()-tick)/((double)CLOCKS_PER_SEC);

  }

printf("\ncp_time = %f\n",cp_time);
printf("mm_time = %f\n",mm_time);
printf("ad_time = %f\n\n",ad_time);


  // return 
  return 0;

}

