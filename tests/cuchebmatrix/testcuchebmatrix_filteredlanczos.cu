#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("../matrices/H2O.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  //string mtxfile("../matrices/G2_circuit.mtx");
  string mtxfile("../matrices/Trefethen_20000.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cucheblanczos
  cucheblanczos ccl;

  // call filtered lanczos
  cuchebmatrix_filteredlanczos(4,0,&ccm,&ccl);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCL
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
