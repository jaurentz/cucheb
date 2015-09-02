#include <cuchebpoly.h>

/* routine for destroying cuchebpoly object */
int cuchebpoly_destroy(cuchebpoly* ccp){

  // free cufft
  cufftDestroy(ccp->cuffthandle);
 
  // free workspace
  cudaFree(ccp->dinput);
  cudaFree(ccp->doutput);
 
  // return 
  return 0;

}

