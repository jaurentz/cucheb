#include <cucheb.h>

/* arnoldi run using cuchebmatrix */
int cucheblanczos_arnoldi(int nsteps, cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, bsize, nblocks, nvecs, stop;
  double scl, one = 1.0, zero = 0.0, mone = -1.0;
  double* bands;
  double* dtemp;
  double* dvecs;
  double* dschurvecs;
  n = ccl->n;
  bsize = ccl->bsize;
  nblocks = ccl->nblocks;
  nvecs = bsize*nblocks;
  stop = ccl->stop;
  bands = ccl->bands;
  dtemp = ccl->dtemp;
  dvecs = ccl->dvecs;
  dschurvecs = ccl->dschurvecs;

  // set niters
  int niters;
  niters = min(nsteps,nblocks-stop);

  // loop through nblocks
  int ind;
  for(int ii=0; ii < niters; ii++){

    // inner loop for bsize blocks
    for(int jj=0; jj < bsize; jj++){

      // set index
      ind = (ii+stop)*bsize + jj;

      // apply matrix
      cuchebmatrix_mv(ccm,&one,&dvecs[ind*n],&zero,&dvecs[(ind+bsize)*n]);

      // orthogonalize
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, (ind+bsize), &one, &dvecs[0], n, 
                  &dvecs[(ind+bsize)*n], 1, &zero,
                  &dschurvecs[ind*(nvecs+bsize)], 1);
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, (ind+bsize), &mone, &dvecs[0], n, 
                  &dschurvecs[ind*(nvecs+bsize)], 1, &one, &dvecs[(ind+bsize)*n], 1);

      // reorthogonalize
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, (ind+bsize), &one, &dvecs[0], n, 
                &dvecs[(ind+bsize)*n], 1, &zero, &dtemp[0], 1);
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, (ind+bsize), &mone, &dvecs[0], n, 
                &dtemp[0], 1, &one, &dvecs[(ind+bsize)*n], 1);
      cublasDaxpy(ccm->cublashandle, (ind+bsize), &one, &dtemp[0], 1, 
                &dschurvecs[ind*(nvecs+bsize)], 1);

      // normalize
      cublasDnrm2(ccm->cublashandle, n, &dvecs[(ind+bsize)*n], 1, &scl);
      cudaMemcpy(&dschurvecs[ind*(nvecs+bsize)+ind+bsize], &scl,
                 sizeof(double), cudaMemcpyHostToDevice);
      scl = 1.0/scl;
      cublasDscal(ccm->cublashandle, n, &scl, &dvecs[(ind+bsize)*n], 1);

    }

  }

  // update stop 
  ccl->stop += niters;

  // copy data to host
  cudaMemcpy(&(ccl->schurvecs)[0],&dschurvecs[0],
             (ccl->stop)*bsize*(nvecs+bsize)*sizeof(double),
             cudaMemcpyDeviceToHost);

  // return  
  return 0;

}

