#include <cuchebmatrix.h>

/* routine for mv multiply on GPU */
int cuchebmatrix_mv(cuchebmatrix* ccm, double* alpha, double* x, double* beta,
                    double* y){

  // cusparseDcsrmv
  cusparseDcsrmv(ccm->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, ccm->m, ccm->n, 
                 ccm->nnz, alpha, ccm->matdescr, ccm->dvals, ccm->drowinds,
                 ccm->dcolinds, x, beta, y);
 
  // return 
  return 0;

}

