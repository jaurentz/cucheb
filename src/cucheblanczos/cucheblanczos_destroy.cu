#include <cucheb.h>

/* routine to free memory in cucheblanczos object */
int cucheblanczos_destroy(cucheblanczos* ccl){

  // free index
  delete[] ccl->index;

  // free diag
  delete[] ccl->diag;

  // free sdiag
  delete[] ccl->sdiag;

  // free schurvecs
  delete[] ccl->schurvecs;

  // free dtemp
  cudaFree(ccl->dtemp);

  // free dvecs
  cudaFree(ccl->dvecs);

  // free dschurvecs
  cudaFree(ccl->dschurvecs);

  // return  
  return 0;

}
