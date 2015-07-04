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

/* routine for standard print */
int cuchebpoly_print(cuchebpoly* ccp){

  // print banner
  printf("\ncuchebpoly:\n");

  // degree
  printf(" degree = %d\n",ccp->degree);
 
  // a and b
  printf(" [a,b] = [%+e,%+e]\n",ccp->a,ccp->b);
  
  // coeffs
  for (int ii=0; ii<ccp->degree+1; ii++) {
    printf(" coeffs[%d] = %+e\n",ii,ccp->coeffs[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

/* routine for long print */
int cuchebpoly_printlong(cuchebpoly* ccp){

  // print banner
  printf("\ncuchebpoly:\n");

  // degree
  printf(" degree = %d\n",ccp->degree);
 
  // a and b
  printf(" [a,b] = [%+e,%+e]\n",ccp->a,ccp->b);
  
  // coeffs
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    printf(" coeffs[%d] = %+e\n",ii,ccp->coeffs[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

/* routine for Chebyshev points */
int cuchebpoints(double a, double b, double* points){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set points
  double alpha = (b-a)/2.0;
  double beta = (b+a)/2.0;
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    points[ii] = alpha*sin(DOUBLE_PI*(2.0*ii-DOUBLE_DEG)/(2.0*DOUBLE_DEG)) +
                  beta;
  }

  // return 
  return 0;

}

/* cuchebcoeffs */
int cuchebcoeffs (double *coeffs){
 
  // allocate workspace
  cufftDoubleReal *input;
  cudaMalloc(&input,2*DOUBLE_DEG*sizeof(cufftDoubleReal));
  cufftDoubleComplex *output;
  cudaMalloc(&output,(DOUBLE_DEG+1)*sizeof(cufftDoubleComplex));
 
  // initialize cufft
  cufftHandle cufftHand;
  cufftPlan1d(&cufftHand, 2*DOUBLE_DEG, CUFFT_D2Z, 1);
 
  // initialize input 
  int deg = DOUBLE_DEG;
  int N = 2*deg;

  cudaMemcpy(&coeffs[deg], &input[0], sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(&coeffs[0], &input[deg], sizeof(double), cudaMemcpyHostToDevice);
  for (int ii=1; ii<deg; ii++) {
    cudaMemcpy(&coeffs[deg-ii], &input[ii], sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&coeffs[deg-ii], &input[N-ii], sizeof(double), cudaMemcpyHostToDevice);
  }
 
  // execute plan
  cufftExecD2Z(cufftHand,input,output);
 
  // extract output
  cudaMemcpy(&output[deg], &coeffs[0], sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&output[0], &coeffs[deg], sizeof(double), cudaMemcpyDeviceToHost);
  for (int ii=1; ii<deg; ii++) {
    cudaMemcpy(&output[deg-ii], &coeffs[ii], sizeof(double), cudaMemcpyDeviceToHost);
  }

  // normalize output
  coeffs[0] = coeffs[0]/(double)(deg)/2.0;
  coeffs[deg] = coeffs[deg]/(double)(deg)/2.0;
  for (int ii=1; ii<deg; ii++) {
    coeffs[ii] = coeffs[ii]/(double)(deg);
  }

  // free cufft
  cufftDestroy(cufftHand);
 
  // free workspace
  cudaFree(input);
  cudaFree(output);
 
  // return success
  return 0;

}
