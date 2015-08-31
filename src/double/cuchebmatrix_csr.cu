#include <cuchebmatrix.h>

/* routine for converting to csr format */
int cuchebmatrix_csr(cuchebmatrix* ccm){

  // loop through row inds
  int cind = 0;
  for(int ii=0; ii<(ccm->nnz); ii++){
    if((ccm->rowinds)[ii] > cind){
      cind += 1;
      (ccm->rowinds)[cind] = ii;
    }
  }
  (ccm->rowinds)[ccm->m] = ccm->nnz;

  // return 
  return 0;

}


