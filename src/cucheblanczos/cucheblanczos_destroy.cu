#include <cucheb.h>
/*
  cucheblanczos_destroy

  This routine frees memory associated with an instance of a cucheblanczos
  object. The following inputs are required:

    ccl - a reference to an instance of a cucheblanczos

*/

/* routine to free memory in cucheblanczos object */
int cucheblanczos_destroy(cucheblanczos* ccl){

  // free index
  delete[] ccl->index;

  // free bands 
  delete[] ccl->bands;

  // free evals
  delete[] ccl->evals;

  // free res
  delete[] ccl->res;

  // free vecs
  delete[] ccl->vecs;

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
