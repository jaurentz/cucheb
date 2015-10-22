#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // cuhebmatrix
  string mtxfile("../matrices/ca2010.mtx");
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // call lanczos for an interval
  cucheblanczos ccl;
  cuchebstats ccstats;
  cuchebmatrix_lanczos(2.0e6, 2.1e6, 1, 1200, 30, &ccm, &ccl, &ccstats);

  // print ccm
  cuchebmatrix_print(&ccm);

  // print ccstats
  cuchebstats_print(&ccstats);

  // print eigenvalues
  for (int ii=0; ii<ccl.nconv; ii++) {
  //for (int ii=0; ii<100; ii++) {
    printf(" %+e, %e\n",ccl.evals[ccl.index[ii]],ccl.res[ccl.index[ii]]);
  }
  printf("\n");

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCB
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
