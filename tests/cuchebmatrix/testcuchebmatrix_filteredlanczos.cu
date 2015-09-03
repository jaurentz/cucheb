#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("./matrices/Si34H36.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);
  cuchebmatrix_print(&ccm);

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
