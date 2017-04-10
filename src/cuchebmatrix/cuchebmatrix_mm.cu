#include <cucheb.h>

/* routine for mm multiply on GPU */
int cuchebmatrix_mm(cuchebmatrix* ccm, int bsize, double* alpha, double* X, 
                    double* beta, double* Y){

  // cusparseDcsrmm
  cusparseDcsrmm(ccm->cusparsehandle, CUSPARSE_OPERATION_NON_TRANSPOSE, ccm->m, bsize, 
                 ccm->n, ccm->nnz, alpha, ccm->matdescr, ccm->dvals, ccm->drowinds,
                 ccm->dcolinds, X, ccm->n, beta, Y, ccm->m);

  // return 
  return 0;

}

