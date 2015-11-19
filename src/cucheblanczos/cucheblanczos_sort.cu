#include <cucheb.h>

/* my interval norm */
double interval_norm(double lb, double ub, double val){

  if (val >= lb && val <= ub) { return val; }
  else { return abs(val) + ub; }

}

/* sort evals in interval */
int cucheblanczos_sort(double lb, double ub, cucheblanczos* ccl){

  // local variables
  int neig;
  int* index;
  double* evals;
  neig = (ccl->nconv);
  index = ccl->index;
  evals = ccl->evals;

  // sort ritz values
  // create a vector of evals and indices
  vector< pair< double , int > > temp;
  for(int ii=0; ii < neig; ii++){
    temp.push_back(make_pair( interval_norm(lb,ub,evals[index[ii]]), index[ii] ));
  }

  // sort vector
  sort(temp.begin(),temp.end());
 
  // update index
  for(int ii=0; ii < neig; ii++){
    index[ii] = temp[ii].second;
  }

  // compute number of converged eigenvalues
  ccl->nconv = 0;
  for(int ii=0; ii < neig; ii++){
    if (evals[index[ii]] > ub || evals[index[ii]] < lb){ break; }
    ccl->nconv += 1;  
  }

  // return  
  return 0;

}



/* sort evals by largest modulus*/
int cucheblanczos_sort(cucheblanczos* ccl){

  // local variables
  int neig;
  int* index;
  double* evals;
  neig = (ccl->nconv);
  index = ccl->index;
  evals = ccl->evals;

  // sort ritz values
  // create a vector of evals and indices
  vector< pair< double , int > > temp;
  for(int ii=0; ii < neig; ii++){
    temp.push_back(make_pair( -evals[index[ii]], index[ii] ));
  }

  // sort vector
  sort(temp.begin(),temp.end());

  // update index
  for(int ii=0; ii < neig; ii++){
    index[ii] = temp[ii].second;
  }

  // return  
  return 0;

}
