#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("../matrices/H2O.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  string mtxfile("../matrices/G2_circuit.mtx");
  //string mtxfile("../matrices/Trefethen_20000.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // compute spectral interval
  cuchebmatrix_specint(&ccm);
  cuchebmatrix_print(&ccm);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
