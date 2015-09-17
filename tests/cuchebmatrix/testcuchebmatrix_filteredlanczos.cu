#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("../matrices/SiH4.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  //string mtxfile("../matrices/H2O.mtx");
  //string mtxfile("../matrices/Si34H36.mtx");
  //string mtxfile("../matrices/Si87H76.mtx");
  //string mtxfile("../matrices/CO.mtx");
  string mtxfile("../matrices/Ga41As41H72.mtx");
  //string mtxfile("../matrices/dielFilterV2real.mtx");
  //string mtxfile("../matrices/G2_circuit.mtx");
  //string mtxfile("../matrices/Trefethen_20000.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cucheblanczos
  cucheblanczos ccl;

  // cuchebstats
  cuchebstats ccstats;

  // call filtered lanczos for a point
  cuchebmatrix_filteredlanczos(10, -1e100, 3, &ccm, &ccl, &ccstats);

  // call filtered lanczos for an interval
  //cuchebmatrix_filteredlanczos(-10.0, -1.15, 3, &ccm, &ccl, &ccstats);

  // print ccm
  cuchebmatrix_print(&ccm);

  // print ccstats
  cuchebstats_print(&ccstats);

  // print eigenvalues
  for (int ii=0; ii<ccstats.num_conv; ii++) {
  //for (int ii=0; ii<ccl.stop*ccl.bsize; ii++) {
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
