#include <cucheb.h>

/* compute ritz values and vectors */
int cuchebblocklanczos_rayleigh(cuchebmatrix* ccm, cuchebblocklanczos* ccb){

  // local variables
  int n, nvecs;
  double* evals;
  double* res;
  double* dtemp;
  double* dvecs;
  n = ccb->n;
  nvecs = (ccb->bsize)*(ccb->nblocks);
  evals = ccb->evals;
  res = ccb->res;
  dtemp = ccm->dtemp;
  dvecs = ccb->dvecs;

  // compute rayleigh quotients and residuals
  double one = 1.0, zero = 0.0;
  double scl, rval;
  for(int ii=0; ii<nvecs; ii++){
 
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
