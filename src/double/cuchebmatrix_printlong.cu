#include <cuchebmatrix.h>

/* routine for long print */
int cuchebmatrix_printlong(cuchebmatrix* ccm){

  // print banner
  printf("\ncuchebmatrix:\n");

  // print matcode
  printf(" matcode = %.4s\n",ccm->matcode);
 
  // print m
  printf(" m = %d\n",ccm->m);
 
  // print n
  printf(" n = %d\n",ccm->n);
 
  // print nnz
  printf(" nnz = %d\n",ccm->nnz);

  // print rowinds, colinds and vals
  for (int ii=0; ii<ccm->nnz; ii++) {
    printf(" rowinds[%d] = %d, colinds[%d] = %d, vals[%d] = %+e\n",
           ii,(ccm->rowinds)[ii],ii,(ccm->colinds)[ii],ii,(ccm->vals)[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

