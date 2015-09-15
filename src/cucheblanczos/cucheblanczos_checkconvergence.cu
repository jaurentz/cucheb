#include <cucheb.h>

/* check convergence of eigenvalues and sort using rho*/
int cucheblanczos_checkconvergence(int* numconv, double rho, cuchebmatrix* ccm, 
                                        cucheblanczos* ccl){

  // local variables
  int nvecs;
  double nrm;
  int* index;
  double* evals;
  double* res;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  nrm = (ccl->n)*max(abs(ccm->a),abs(ccm->b));
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;

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

/* helper routine for sorting */
double mynrm(double lb, double ub, double num){

  if ( num >= lb && num <= ub ) { return num; }
  else { return (ub + abs(num-(ub-lb)/2.0)); }

}

/* check convergence of eigenvalues and sort smallest to largest */
int cucheblanczos_checkconvergence(int* numconv, double lb, double ub, cuchebmatrix* ccm, 
                                        cucheblanczos* ccl){

  // local variables
  int nvecs;
  double nrm;
  int* index;
  double* evals;
  double* res;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  nrm = (ccl->n)*max(abs(ccm->a),abs(ccm->b));
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;

  // compute number of converged eigenvalues
  *numconv = 0;
  for(int ii=0; ii < nvecs; ii++){

    if (res[index[ii]] >= DOUBLE_TOL*nrm) {
      *numconv = ii;
      break;
    }

  }

  // create a vector of evals and indices
  vector< pair< double , int > > temp;
  for(int ii=0; ii < *numconv; ii++){
    temp.push_back(make_pair( mynrm(lb,ub,evals[index[ii]]) , index[ii] ));
  }

  // sort vector
  sort(temp.begin(),temp.end());

  // update index
  for(int ii=0; ii < *numconv; ii++){
    index[ii] = temp[ii].second;
  }

  // update numconv
  int cnt = 0;
  for(int ii=0; ii < *numconv; ii++){
    if (evals[index[ii]] >= lb && evals[index[ii]] <= ub) { cnt += 1; }
    else { break; }
  }
  
  // set numconv
  *numconv = cnt;

  // return  
  return 0;

}
