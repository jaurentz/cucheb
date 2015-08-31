#include <cucheblanczos.h>

/* routine to free memory in cucheblanczos object */
int cucheblanczos_destroy(cucheblanczos* ccl){

  // free diag
  delete[] ccl->diag;

  // free sdiag
  delete[] ccl->sdiag;

  // free schurvecs
  delete[] ccl->schurvecs;

  // destroy cublas handle
  cublasDestroy(ccl->handle);
 
  // free dvecs
  cudaFree(ccl->dvecs);

  // free dschurvecs
  cudaFree(ccl->dschurvecs);

  // return  
  return 0;

}
