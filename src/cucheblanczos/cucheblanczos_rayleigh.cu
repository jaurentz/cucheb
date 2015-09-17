#include <cucheb.h>

/* compute ritz values and vectors */
int cucheblanczos_rayleigh(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, bsize, nvecs, stop;
  int* index;
  double* evals;
  double* res;
  double* vecs;
  double* schurvecs;
  double* dv1;
  double* dv2;
  double* dschurvecs;
  double* dvecs;
  n = ccl->n;
  bsize = ccl->bsize;
  stop = ccl->stop;
  nvecs = bsize*(ccl->nblocks);
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;
  schurvecs = ccl->schurvecs;
  vecs = ccl->vecs;
  dv1 = &(ccm->dtemp)[0];
  dv2 = &(ccm->dtemp)[n];
  dschurvecs = ccl->dschurvecs;
  dvecs = ccl->dvecs;

  // copy schurvecs into dschur
  for(int ii=0; ii<stop*bsize; ii++){
    cudaMemcpy(&dschurvecs[ii*(nvecs+bsize)],&schurvecs[index[ii]*(nvecs+bsize)],
               (nvecs+bsize)*sizeof(double),cudaMemcpyHostToDevice);
  }

  // compute rayleigh quotients and residuals
  double one = 1.0, zero = 0.0;
  double scl, rval;
  for(int ii=0; ii<stop*bsize; ii++){
 
    // compute ritz vector
    cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, stop*bsize, &one, dvecs, 
                n, &dschurvecs[ii*(nvecs+bsize)], 1, &zero, dv1, 1);

    // copy ritz vector to cpu
    cudaMemcpy(&vecs[ii*n],dv1,n*sizeof(double),cudaMemcpyDeviceToHost);

    // apply operator
    cuchebmatrix_mv(ccm,&one,dv1,&zero,dv2);

    // compute rayleigh quotient
    cublasDnrm2(ccm->cublashandle,n,dv1,1,&scl);
    cublasDdot(ccm->cublashandle,n,dv1,1,dv2,1,&evals[ii]);
    evals[ii] = evals[ii]/scl/scl;

    // compute residual vector
    rval = -evals[ii];
    cublasDaxpy(ccm->cublashandle,n,&rval,dv1,1,dv2,1);

    // compute norm of residual
    cublasDnrm2(ccm->cublashandle,n,dv2,1,&res[ii]);
    res[ii] = res[ii]/scl;

    // reset index 
    index[ii] = ii;
  
  }

  // return  
  return 0;

}
