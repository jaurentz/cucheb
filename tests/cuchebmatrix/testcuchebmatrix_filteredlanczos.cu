#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("./matrices/Sandi_authors.mtx");
  //string mtxfile("./matrices/Trefethen_20.mtx");
  //string mtxfile("./matrices/Stranke94.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);
  cuchebmatrix_print(&ccm);

  // cucheblanczos
  cucheblanczos ccl;

  // call filtered lanczos
  cuchebmatrix_filteredlanczos(4,0,&ccm,&ccl);

  // print eigenvalues
  for(int ii=0; ii < ccl.nvecs; ii++){
    printf(" diag[%d] = %+e, sdiag[%d] = %+e\n",
           ii,ccl.diag[ii],ii,ccl.sdiag[ii]);
  }
  printf("\n");

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCL
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
