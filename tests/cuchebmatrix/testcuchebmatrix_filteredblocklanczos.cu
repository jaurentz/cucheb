#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("../matrices/H2O.mtx");
  //string mtxfile("../matrices/Ga41As41H72.mtx");
  //string mtxfile("../matrices/Si87H76.mtx");
  string mtxfile("../matrices/dielFilterV2real.mtx");
  //string mtxfile("../matrices/CO.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  //string mtxfile("../matrices/G2_circuit.mtx");
  //string mtxfile("../matrices/Trefethen_20000.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cuchebblocklanczos
  cuchebblocklanczos ccb;

  // call filtered lanczos
  cuchebmatrix_filteredblocklanczos(4,-1e100,3,&ccm,&ccb);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCB
  cuchebblocklanczos_destroy(&ccb);

  // return 
  return 0;

}
