#include <cucheb.h>

/* compute ritz values and vectors */
int cucheblanczos_rayleigh(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, bsize, nvecs, stop;
  int* index;
  double* evals;
  double* res;
  double* bands;
  double* schurvecs;
  double* dtemp;
  double* dschurvecs;
  double* dvecs;
  n = ccl->n;
  bsize = ccl->bsize;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  stop = ccl->stop;
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;
  bands = ccl->bands;
  schurvecs = ccl->schurvecs;
  dtemp = ccm->dtemp;
  dschurvecs = ccl->dschurvecs;
  dvecs = ccl->dvecs;

  // reset index
  for(int ii=0; ii < stop*bsize; ii++){
    index[ii] = ii;
  }

  // fill bands
  for(int jj=0; jj<stop*bsize; jj++) {
    for(int ii=0; ii<bsize; ii++){
      bands[jj*(bsize+1)+ii] = schurvecs[jj*(nvecs+bsize)+jj+ii];
    }
  }

  // initialize residuals
  for(int ii=0; ii < stop*bsize; ii++){
    res[ii] = bands[(bsize+1)*stop*bsize-1];
  }

  // initialize schurvectors
  for(int jj=0; jj<stop*bsize; jj++) {
    for(int ii=0; ii<nvecs+bsize; ii++){
      if (ii == jj){ schurvecs[jj*(nvecs+bsize)+ii] = 1.0; }
      else{ schurvecs[jj*(nvecs+bsize)+ii] = 0.0; }
    }
  }

  // call bandsymqr
  cuchebutils_bandsymqr(stop*bsize, bsize+1, bands, bsize+1,
                evals, schurvecs, nvecs+bsize);

  // update schurvecs and dschurvecs
  for(int ii=0; ii < stop*bsize; ii++){
    
res[ii] = res[ii]*abs(schurvecs[ii*(nvecs+bsize)+stop*bsize-1]);
printf(" %+e, %e\n",evals[ii],res[ii]);

    // copy schurvecs into dtemp
    cudaMemcpy(&dtemp[0],&schurvecs[ii*(nvecs+bsize)],
               (nvecs+bsize)*sizeof(double),cudaMemcpyHostToDevice);

    // update dschurvecs into schurvecs
    cudaMemcpy(&schurvecs[ii*(nvecs+bsize)],&dschurvecs[ii*(nvecs+bsize)],
               (nvecs+bsize)*sizeof(double),cudaMemcpyDeviceToHost);

    // copy dtemp into dschurvecs
    cudaMemcpy(&dschurvecs[ii*(nvecs+bsize)],&dtemp[0],
               (nvecs+bsize)*sizeof(double),cudaMemcpyDeviceToDevice);

  }

  // update dvecs
  double one = 1.0, zero = 0.0;
  cublasDgemm(ccm->cublashandle, CUBLAS_OP_N, CUBLAS_OP_N, n, stop*bsize,
              stop*bsize, &one, &dvecs[0], n, &dschurvecs[0], nvecs+bsize, 
              &zero, &dvecs[0], n);

  // compute rayleigh quotients and residuals
  double scl, rval;
  for(int ii=0; ii<stop*bsize; ii++){
 
    // apply operator
    cuchebmatrix_mv(ccm,&one,&dvecs[ii*n],&zero,dtemp);

    // compute rayleigh quotient
    cublasDnrm2(ccm->cublashandle,n,&dvecs[ii*n],1,&scl);
    cublasDdot(ccm->cublashandle,n,dtemp,1,&dvecs[ii*n],1,&evals[ii]);
    evals[ii] = evals[ii]/scl/scl;

    // compute residual vector
    rval = -evals[ii];
    cublasDaxpy(ccm->cublashandle,n,&rval,&dvecs[ii*n],1,dtemp,1);

    // compute norm of residual
    cublasDnrm2(ccm->cublashandle,n,dtemp,1,&res[ii]);
    res[ii] = res[ii]/scl;
  
  }

  // return  
  return 0;

}
