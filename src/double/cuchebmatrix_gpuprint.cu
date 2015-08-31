#include <cuchebmatrix.h>

/* routine for long print */
int cuchebmatrix_gpuprint(cuchebmatrix* ccm){

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
  printf("\n");

  // print rowinds
  int ind;
  for (int ii=0; ii<(ccm->m)+1; ii++) {
    cudaMemcpy(&ind,&(ccm->drowinds)[ii],sizeof(int),cudaMemcpyDeviceToHost);
    printf(" rowinds[%d] = %d\n",ii,ind);
  }
  printf("\n");

  // print colinds and vals
  double val;
  for (int ii=0; ii<ccm->nnz; ii++) {
    cudaMemcpy(&ind,&(ccm->dcolinds)[ii],sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(&val,&(ccm->dvals)[ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" colinds[%d] = %d, vals[%d] = %+e\n",ii,ind,ii,val);
  }
  printf("\n");

  // return 
  return 0;

}

