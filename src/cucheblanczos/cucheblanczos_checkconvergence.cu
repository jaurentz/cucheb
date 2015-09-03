#include <cucheb.h>

/* check convergence of eigenvalues */
int cucheblanczos_checkconvergence(int* numconv, double rho, cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, nvecs;
  double nrm;
  int* index;
  double* diag;
  double* sdiag;
  n = ccl->n;
  nvecs = ccl->nvecs;
  nrm = n*max(abs(ccm->a),abs(ccm->b));
  index = ccl->index;
  diag = ccl->diag;
  sdiag = ccl->sdiag;

  // compute number of converged eigenvalues
  *numconv = 0;
  for(int ii=0; ii < nvecs; ii++){

    if (sdiag[index[ii]] >= DOUBLE_TOL*nrm) {
      *numconv = ii;
      break;
    }

  }

  // sort converged ones based on distance from rho
  // create a vector of evals and indices
  vector< pair< double , int > > evals;
  for(int ii=0; ii < *numconv; ii++){
    evals.push_back(make_pair( abs(diag[ii]-rho) , index[ii] ));
  }

  // sort vector
  sort(evals.begin(),evals.end());

  // update index
  for(int ii=0; ii < *numconv; ii++){
    index[ii] = evals[ii].second;
  }

  // return  
  return 0;

}
