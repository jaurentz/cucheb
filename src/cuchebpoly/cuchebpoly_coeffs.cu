#include <cucheb.h>

/* cuchebpoly_coeffs */
int cuchebpoly_coeffs (cuchebpoly* ccp){
 
  // initialize input 
  int deg = DOUBLE_DEG;
  cudaMemcpy(&(ccp->dinput)[0], &(ccp->points)[0], 2*deg*sizeof(double), cudaMemcpyHostToDevice);

  // execute plan
  cufftExecD2Z(ccp->cuffthandle,ccp->dinput,ccp->doutput);
 
  // extract output
  for (int ii=0; ii < deg+1; ii++) {
    cudaMemcpy(&(ccp->coeffs)[ii], &(ccp->doutput)[ii], sizeof(double), cudaMemcpyDeviceToHost);
  }

  // normalize output
  (ccp->coeffs)[0] = (ccp->coeffs)[0]/(double)(deg)/2.0;
  (ccp->coeffs)[deg] = (ccp->coeffs)[deg]/(double)(deg)/2.0;
  for (int ii=1; ii<deg; ii++) {
    (ccp->coeffs)[ii] = (ccp->coeffs)[ii]/(double)(deg);
  }

  // return success
  return 0;

}

