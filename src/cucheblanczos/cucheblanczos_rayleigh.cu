#include <cucheb.h>

/* compute ritz values and vectors */
int cucheblanczos_rayleigh(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, nvecs;
  double* diag;
  double* sdiag;
  double* dtemp;
  double* dvecs;
  n = ccl->n;
  nvecs = ccl->nvecs;
  diag = ccl->diag;
  sdiag = ccl->sdiag;
  dtemp = ccm->dtemp;
  dvecs = ccl->dvecs;

  // compute rayleigh quotients and residuals
  double one = 1.0, zero = 0.0;
  double scl, rval;
  for(int ii=0; ii<nvecs; ii++){
 
    // apply operator
    cuchebmatrix_mv(ccm,&one,&dvecs[ii*n],&zero,dtemp);

    // compute rayleigh quotient
    cublasDnrm2(ccm->cublashandle,n,&dvecs[ii*n],1,&scl);
    cublasDdot(ccm->cublashandle,n,dtemp,1,&dvecs[ii*n],1,&diag[ii]);
    diag[ii] = diag[ii]/scl/scl;

    // compute residual vector
    rval = -diag[ii];
    cublasDaxpy(ccm->cublashandle,n,&rval,&dvecs[ii*n],1,dtemp,1);

    // compute norm of residual
    cublasDnrm2(ccm->cublashandle,n,dtemp,1,&sdiag[ii]);
    sdiag[ii] = sdiag[ii]/scl;

  }


  // return  
  return 0;

}
