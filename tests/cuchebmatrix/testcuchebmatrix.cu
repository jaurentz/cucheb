#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("../matrices/Stranke94.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;

  // initialize CCM
  cuchebmatrix_init(mtxfile, &ccm);

  // print CCM
  cuchebmatrix_print(&ccm);

  // printlong CCM
  cuchebmatrix_printlong(&ccm);

  // printlong CCM
  cuchebmatrix_gpuprint(&ccm);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
