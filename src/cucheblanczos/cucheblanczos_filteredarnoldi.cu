#include <cucheb.h>

/* arnoldi run using cuchebmatrix */
int cucheblanczos_filteredarnoldi(cuchebmatrix* ccm, cuchebpoly* ccp,
                                  cucheblanczos* ccl, cuchebstats* ccstats){

  // local variables
  int n, bsize, nblocks, nvecs, odepth;
  double scl, one = 1.0, zero = 0.0, mone = -1.0;
  double* bands;
  double* dtemp;
  double* dvecs;
  double* dschurvecs;
  n = ccl->n;
  bsize = ccl->bsize;
  nblocks = ccl->nblocks;
  nvecs = bsize*nblocks;
  bands = ccl->bands;
  dtemp = ccl->dtemp;
  dvecs = ccl->dvecs;
  dschurvecs = ccl->dschurvecs;

  // loop through nblocks
  int ind, start;
  for(int ii=0; ii < nblocks; ii++){

    // inner loop for bsize blocks
    for(int jj=0; jj < bsize; jj++){

      // set index
      ind = ii*bsize + jj;

      // apply matrix
      cuchebmatrix_polymv(ccm,ccp,&dvecs[ind*n],&dvecs[(ind+bsize)*n]);

      // num_matvecs
      ccstats->num_matvecs += (ccp->degree);

      // compute orthogonalization depth
      odepth = min((MAX_ORTH_DEPTH)*bsize+jj,ind+bsize);
      start = ind + bsize - odepth;

      // orthogonalize
//      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, (ind+bsize), &one, &dvecs[0], n, 
//                  &dvecs[(ind+bsize)*n], 1, &zero, &dschurvecs[ind*(nvecs+bsize)], 1);
//      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, (ind+bsize), &mone, &dvecs[0], n, 
//                  &dschurvecs[ind*(nvecs+bsize)], 1, &one, &dvecs[(ind+bsize)*n], 1);

      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, odepth, &one, &dvecs[start*n], 
                  n, &dvecs[(ind+bsize)*n], 1, &zero, &dschurvecs[ind*(nvecs+bsize)+start], 1);
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, odepth, &mone, &dvecs[start*n], 
                  n, &dschurvecs[ind*(nvecs+bsize)+start], 1, &one, &dvecs[(ind+bsize)*n], 1);

      // num_innerprods 
      ccstats->num_innerprods += odepth;

      // reorthogonalize
//      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, (ind+bsize), &one, &dvecs[0], n, 
//                &dvecs[(ind+bsize)*n], 1, &zero, &dtemp[0], 1);
//      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, (ind+bsize), &mone, &dvecs[0], n, 
//                &dtemp[0], 1, &one, &dvecs[(ind+bsize)*n], 1);
//      cublasDaxpy(ccm->cublashandle, (ind+bsize), &one, &dtemp[0], 1, 
//                &dschurvecs[ind*(nvecs+bsize)], 1);

      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, odepth, &one, &dvecs[start*n], 
                  n, &dvecs[(ind+bsize)*n], 1, &zero, &dtemp[0], 1);
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, odepth, &mone, &dvecs[start*n], 
                  n, &dtemp[0], 1, &one, &dvecs[(ind+bsize)*n], 1);
      cublasDaxpy(ccm->cublashandle, odepth, &one, &dtemp[0], 1, 
                  &dschurvecs[ind*(nvecs+bsize)+start], 1);

      // num_innerprods 
      ccstats->num_innerprods += odepth;

      // normalize
      cublasDnrm2(ccm->cublashandle, n, &dvecs[(ind+bsize)*n], 1,
                  &bands[(ind+1)*(bsize+1)-1]);
      scl = 1.0/bands[(ind+1)*(bsize+1)-1];
      cublasDscal(ccm->cublashandle, n, &scl, &dvecs[(ind+bsize)*n], 1);

    }

  }

  // copy data to host
  cudaMemcpy(&(ccl->schurvecs)[0],&dschurvecs[0],nvecs*(nvecs+bsize)*sizeof(double),
             cudaMemcpyDeviceToHost);

  // return  
  return 0;

}

