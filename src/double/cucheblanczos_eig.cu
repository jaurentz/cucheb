#include <cucheblanczos.h>

/* compute ritz values and vectors */
int cucheblanczos_eig(cucheblanczos* ccl){

  // local variables
  int nvecs;
  double* diag;
  double* sdiag;
  double* schurvecs;
  nvecs = ccl->nvecs;
  diag = ccl->diag;
  sdiag = ccl->sdiag;
  schurvecs = ccl->schurvecs;

  // fill diag
  for(int ii=0; ii<nvecs; ii++){
    diag[ii] = schurvecs[ii*nvecs + ii];
  }

  // call lapack
  LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'i', nvecs, &diag[0], &sdiag[0],
                 &schurvecs[0], nvecs);

  // return  
  return 0;

}
