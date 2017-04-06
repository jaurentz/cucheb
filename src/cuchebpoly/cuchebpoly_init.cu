#include <cucheb.h>


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


/* routine for basic initialization */
int cuchebpoly_init(cuchebpoly* ccp){

  // set degree
  ccp->degree = 0;
  
  // set a and b
  ccp->a = -1.0;
  ccp->b = 1.0;
  
  // initialize cufftHandle
  cufftPlan1d(&(ccp->cuffthandle), 2*MAX_DOUBLE_DEG, CUFFT_D2Z, 1);

size_t freeMem, totalMem;
cudaMemGetInfo(&freeMem, &totalMem);
printf("cuchebpoly_init\n");
printf("Free = %ld, Total = %ld\n", freeMem, totalMem);

  // allocate workspace
//  if(cudaMalloc(&(ccp->dinput),2*MAX_DOUBLE_DEG*sizeof(cufftDoubleReal)) != 0) {
//    printf("Device memory allocation failed: dinput\n");
gpuErrchk(cudaMalloc(&(ccp->dinput),2*MAX_DOUBLE_DEG*sizeof(cufftDoubleReal)));
//    exit(1);
//  }
  if(cudaMalloc(&(ccp->doutput),(MAX_DOUBLE_DEG+1)*sizeof(cufftDoubleComplex)) != 0) {
    printf("Device memory allocation failed: doutput\n");
    exit(1);
  }

cudaMemGetInfo(&freeMem, &totalMem);
printf("Free = %ld, Total = %ld\n", freeMem, totalMem);
 
  // return 
  return 0;

}

