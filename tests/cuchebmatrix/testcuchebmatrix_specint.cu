#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("../matrices/SiH4.mtx");

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
