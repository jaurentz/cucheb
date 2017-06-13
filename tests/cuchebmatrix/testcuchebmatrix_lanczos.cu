#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(0);

  // cuhebmatrix
  string mtxfile("../matrices/SiH4.mtx");
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // call lanczos for an interval
  cucheblanczos ccl;
  cuchebstats ccstats;
  cuchebmatrix_lanczos(15.0, 16.0, 1, 3000, 3000, &ccm, &ccl, &ccstats);

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
