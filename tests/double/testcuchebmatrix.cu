#include <cuchebmatrix.h>

/* driver */
int main(){

  // cuhebmatrix
  cuchebmatrix ccm;

  // initialize CCM
  cuchebmatrix_init("./matrices/Trefethen_20.mtx",&ccm);

  // print CCM
  cuchebmatrix_print(&ccm);

  // printlong CCM
  cuchebmatrix_printlong(&ccm);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
