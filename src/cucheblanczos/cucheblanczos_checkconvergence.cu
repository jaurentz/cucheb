#include <cucheb.h>

/* check residuals for convergence */
int cucheblanczos_checkconvergence(cucheblanczos* ccl){

  // local variables
  int neig;
  int* index;
  double* evals;
  double* res;
  neig = (ccl->nconv);
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;

  // compute norm of projected matrix
  double nrm = 0.0;
  for(int ii=0; ii < (ccl->bsize)*(ccl->stop); ii++){
    if (abs(evals[index[ii]]) > nrm){ nrm = abs(evals[index[ii]]); }
  }

  // compute number of converged eigenvalues
  nrm = nrm*DOUBLE_TOL;
  ccl->nconv = 0;
  for(int ii=0; ii < neig; ii++){
    if (res[index[ii]] < nrm){ ccl->nconv += 1; }
  }

  // return  
  return 0;

}
