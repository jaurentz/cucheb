#include <cucheb.h>

/* helper routine for sorting */
double mynrm(double lb, double ub, double num){

  if ( num >= lb && num <= ub ) { return num; }
  else { return (ub + abs(num-(ub-lb)/2.0)); }

}

/* check convergence of eigenvalues and sort smallest to largest */
int cucheblanczos_checkconvergence(int* numconv, double lb, double ub, 
                                   cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int nvecs;
  double nrm;
  int* index;
  double* evals;
  double* res;
  nvecs = (ccl->bsize)*(ccl->stop);
  nrm = sqrt(1.0*(ccl->n))*max(abs(ccm->a),abs(ccm->b));
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;

  // sort rayleigh quotients
  // create a vector of evals and indices
  vector< pair< double , int > > temp;
  for(int ii=0; ii < nvecs; ii++){
    temp.push_back(make_pair( mynrm(lb,ub,evals[index[ii]]) , index[ii] ));
  }

  // sort vector
  sort(temp.begin(),temp.end());

  // update index
  for(int ii=0; ii < nvecs; ii++){
    index[ii] = temp[ii].second;
  }

  // compute number of evals in [lb,ub]
  int numint = 0;
  for (int ii=0; ii < nvecs; ii++){
    if (evals[index[ii]] < lb || evals[index[ii]] > ub) { break; }
    else { numint = ii+1; }
  }

  // compute number of converged eigenvalues
  *numconv = 0;
  for(int ii=0; ii < numint; ii++){
    if (res[index[ii]] >= DOUBLE_TOL*nrm){ break; }
    else { *numconv = ii+1; }
  }

  // update numconv
  if (*numconv < numint) { *numconv = 0; }

  // return  
  return 0;

}
