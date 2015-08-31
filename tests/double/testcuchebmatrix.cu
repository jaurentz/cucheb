#include <cuchebmatrix.h>

/* driver */
int main(){

  // input file
  string mtxfile("./matrices/Trefethen_20.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;

  // initialize CCM
  cuchebmatrix_init(mtxfile, &ccm);

  // print CCM
  cuchebmatrix_print(&ccm);

  // printlong CCM
  cuchebmatrix_printlong(&ccm);

  // sort entries of CCM
  cuchebmatrix_sort(&ccm);

  // printlong CCM
  cuchebmatrix_printlong(&ccm);

  // convert CCM to csr format
  cuchebmatrix_csr(&ccm);

  // printlong CCM
  cuchebmatrix_printlong(&ccm);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
