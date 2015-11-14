#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // cuhebmatrix
  //string mtxfile("../matrices/ca2010.mtx");
  //string mtxfile("../matrices/mn2010.mtx");
  //string mtxfile("../matrices/Ge99H100.mtx");
  string mtxfile("../matrices/caidaRouterLevel.mtx");
  //string mtxfile("../matrices/rgg_n_2_20_s0.mtx");
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  double a, b, per;
  a = -108;
  b = 109;
  per = .42;

  // call filtered lanczos for an interval
  cucheblanczos ccl;
  cuchebstats ccstats;
  //cuchebmatrix_filteredlanczos(4.0e6, 4.1e6, 1, &ccm, &ccl, &ccstats);
  //cuchebmatrix_expertlanczos(3.0e6, 3.5e6, -1, 1, 1200, 30, &ccm, &ccl, &ccstats);
  cuchebmatrix_expertlanczos(b-per*abs(b-a), 1.1*b, -1, 1, 1200, 30, &ccm, &ccl, &ccstats);

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
