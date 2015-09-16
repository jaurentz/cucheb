#include <cucheb.h>

/* reset ccl so that more Arnoldi vectors can be added */
int cucheblanczos_reset(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, bsize, nvecs, stop;
  double* schurvecs;
  double* dschurvecs;
  double* dvecs;
  n = ccl->n;
  bsize = ccl->bsize;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  stop = ccl->stop;
  schurvecs = ccl->schurvecs;
  dschurvecs = ccl->dschurvecs;
  dvecs = ccl->dvecs;

  // update dvecs
  double one = 1.0, zero = 0.0;
  cublasDgemm(ccm->cublashandle, CUBLAS_OP_N, CUBLAS_OP_T, n, stop*bsize,
              stop*bsize, &one, &dvecs[0], n, &dschurvecs[0], nvecs+bsize, 
              &zero, &dvecs[0], n);

  // copy schurvecs into dschurvecs
  cudaMemcpy(&dschurvecs[0],&schurvecs[0],
             stop*bsize*(nvecs+bsize)*sizeof(double),cudaMemcpyHostToDevice);

  // return  
  return 0;

}
