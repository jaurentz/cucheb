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
    else { ccl->nconv = ii+1; }
  }

  // return  
  return 0;

}



/* sort evals around point*/
int cucheblanczos_sort(double rho, cucheblanczos* ccl){

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
    temp.push_back(make_pair( abs(rho-evals[index[ii]]), index[ii] ));
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
