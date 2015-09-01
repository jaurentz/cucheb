#include <cucheblanczos.h>

/* driver */
int main(){

  // input file
  //string mtxfile("./matrices/Sandi_authors.mtx");
  //string mtxfile("./matrices/Trefethen_20.mtx");
  string mtxfile("./matrices/Stranke94.mtx");

  // cuchebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cucheblanczos
  cucheblanczos ccl;
  cucheblanczos_init(&ccm, &ccl);

  // print CCL
  cucheblanczos_print(&ccl);

  // set starting vector
  cucheblanczos_startvec(&ccl);

  // do arnoldi run
  cucheblanczos_arnoldi(&ccm,&ccl);

  // print starting vector
  double val;
  for(int jj=0; jj < ccl.nvecs+1; jj++){
  for(int ii=0; ii < ccl.n; ii++){
    cudaMemcpy(&val,&(ccl.dvecs)[jj*ccl.n + ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dvecs[%d] = %+e\n", jj*ccl.n+ii, val);
  }
  printf("\n");
  }
  printf("\n");

  // print schurvecs
  for(int jj=0; jj < ccl.nvecs; jj++){
  for(int ii=0; ii < ccl.nvecs; ii++){
    printf(" schurvecs[%d] = %+e\n", jj*ccl.nvecs+ii, ccl.schurvecs[jj*ccl.nvecs+ii]);
  }
  printf("\n");
  }
  printf("\n");

  // print sdiag
  for(int ii=0; ii < ccl.nvecs; ii++){
    printf(" sdiag[%d] = %+e\n", ii, ccl.sdiag[ii]);
  }
  printf("\n");

  // compute ritz values
  cucheblanczos_eig(&ccl);

  // print diag and sdiag
  for(int ii=0; ii < ccl.nvecs; ii++){
    printf(" diag[%d] = %+e, sdiag[%d] = %+e\n", ii, ccl.diag[ii], ii, ccl.sdiag[ii]);
  }
  printf("\n");

  // print schurvecs
  for(int jj=0; jj < ccl.nvecs; jj++){
  for(int ii=0; ii < ccl.nvecs; ii++){
    printf(" schurvecs[%d] = %+e\n", jj*ccl.nvecs+ii, ccl.schurvecs[jj*ccl.nvecs+ii]);
  }
  printf("\n");
  }
  printf("\n");

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCL
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
