#include <cucheb.h>

/* routine to free memory in cuchebblocklanczos object */
int cuchebblocklanczos_destroy(cuchebblocklanczos* ccb){

  // free index
  delete[] ccb->index;

  // free bands 
  delete[] ccb->bands;

  // free evals
  delete[] ccb->evals;

  // free res
  delete[] ccb->res;

  // free schurvecs
  delete[] ccb->schurvecs;

  // free dtemp
  cudaFree(ccb->dtemp);

  // free dvecs
  cudaFree(ccb->dvecs);

  // free dschurvecs
  cudaFree(ccb->dschurvecs);

  // return  
  return 0;

}
