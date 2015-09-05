#include <cucheb.h>

/* check convergence of eigenvalues */
int cuchebblocklanczos_checkconvergence(int* numconv, double rho, cuchebmatrix* ccm, 
                                        cuchebblocklanczos* ccb){

  // local variables
  int n, nvecs;
  double nrm;
  int* index;
  double* evals;
  double* res;
  n = ccb->n;
  nvecs = (ccb->bsize)*(ccb->nblocks);
  nrm = n*max(abs(ccm->a),abs(ccm->b));
  index = ccb->index;
  evals = ccb->evals;
  res = ccb->res;

  // compute number of converged eigenvalues
  *numconv = 0;
  for(int ii=0; ii < nvecs; ii++){

    if (res[index[ii]] >= DOUBLE_TOL*nrm) {
      *numconv = ii;
      break;
    }

  }

  // sort converged ones based on distance from rho
  // create a vector of evals and indices
  vector< pair< double , int > > temp;
  for(int ii=0; ii < *numconv; ii++){
    temp.push_back(make_pair( abs(evals[ii]-rho) , index[ii] ));
  }

  // sort vector
  sort(temp.begin(),temp.end());

  // update index
  for(int ii=0; ii < *numconv; ii++){
    index[ii] = temp[ii].second;
  }

  // return  
  return 0;

}
