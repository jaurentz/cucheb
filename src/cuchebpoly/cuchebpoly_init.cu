#include <cucheb.h>

/* routine for basic initialization */
int cuchebpoly_init(cuchebpoly* ccp){

  // set degree
  ccp->degree = 0;
  
  // set a and b
  ccp->a = -1.0;
  ccp->b = 1.0;
  
  // initialize cufftHandle
  cufftPlan1d(&(ccp->cuffthandle), 2*MAX_DOUBLE_DEG, CUFFT_D2Z, 1);

  // allocate workspace
  if(cudaMalloc(&(ccp->dinput),2*MAX_DOUBLE_DEG*sizeof(cufftDoubleReal)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccp->doutput),(MAX_DOUBLE_DEG+1)*sizeof(cufftDoubleComplex)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
 
  // return 
  return 0;

}

