#include <cuchebpoly.h>

/* cuchebpoly_coeffs */
int cuchebpoly_coeffs (double *coeffs){
 
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

  cudaMemcpy(&input[0], &coeffs[deg], sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(&input[deg], &coeffs[0], sizeof(double), cudaMemcpyHostToDevice);
  for (int ii=1; ii<deg; ii++) {
    cudaMemcpy(&input[ii], &coeffs[deg-ii], sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&input[N-ii], &coeffs[deg-ii], sizeof(double), cudaMemcpyHostToDevice);
  }

  // execute plan
  cufftExecD2Z(cufftHand,input,output);
 
  // extract output
  cudaMemcpy(&coeffs[0], &output[0], sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&coeffs[deg], &output[deg], sizeof(double), cudaMemcpyDeviceToHost);
  for (int ii=1; ii<deg; ii++) {
    cudaMemcpy(&coeffs[ii], &output[ii], sizeof(double), cudaMemcpyDeviceToHost);
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

