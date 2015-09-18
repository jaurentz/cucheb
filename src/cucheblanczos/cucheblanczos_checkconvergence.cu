#include <cucheb.h>

/* check convergence of eigenvalues */
int cucheblanczos_checkconvergence(int* numconv, cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int nvecs;
  double nrm;
  int* index;
  double* res;
  nvecs = (ccl->bsize)*(ccl->stop);
  nrm = max(abs(ccm->a),abs(ccm->b));
  index = ccl->index;
  res = ccl->res;

  // compute number of converged eigenvalues
  *numconv = 0;
  for(int ii=0; ii < nvecs; ii++){
    if (res[index[ii]] >= DOUBLE_TOL*nrm){ break; }
    else { *numconv = ii+1; }
  }

  // return  
  return 0;

}

