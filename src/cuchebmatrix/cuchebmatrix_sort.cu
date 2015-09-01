#include <cuchebmatrix.h>

/* routine for sorting entries */
int cuchebmatrix_sort(cuchebmatrix* ccm){

  // create a vector of pairs of pairs
  vector< pair< pair<int,int> , double > > mat;
  for(int ii=0; ii<(ccm->nnz); ii++){
    mat.push_back(make_pair( make_pair((ccm->rowinds)[ii],(ccm->colinds)[ii]) , (ccm->vals)[ii] ));
  }

  // sort vector
  sort(mat.begin(),mat.end());

  // update ccm
  for(int ii=0; ii < mat.size(); ii++){
    (ccm->rowinds)[ii] = mat[ii].first.first;
    (ccm->colinds)[ii] = mat[ii].first.second;
    (ccm->vals)[ii] = mat[ii].second;
  }

  // return 
  return 0;

}

