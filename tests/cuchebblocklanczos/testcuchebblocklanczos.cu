#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("./matrices/Sandi_authors.mtx");
  //string mtxfile("../matrices/Trefethen_20.mtx");
  string mtxfile("../matrices/Stranke94.mtx");

  // cuchebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cuchebblocklanczos
  cuchebblocklanczos ccb;
  cuchebblocklanczos_init(1, 10, &ccm, &ccb);

  // print CCB
  cuchebblocklanczos_print(&ccb);

  // set starting vector
  cuchebblocklanczos_startvecs(&ccb);

  // do arnoldi run
  cuchebblocklanczos_arnoldi(&ccm,&ccb);

  // print arnoldi vectors
  double val;
  int nvecs;
  nvecs = (ccb.bsize)*(ccb.nblocks);
  for(int jj=0; jj < nvecs+ccb.bsize; jj++){
  for(int ii=0; ii < ccb.n; ii++){
    cudaMemcpy(&val,&(ccb.dvecs)[jj*ccb.n + ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dvecs[%d] = %+e\n", jj*ccb.n+ii, val);
  }
  printf("\n");
  }
  printf("\n");

  // compute ritz values
  cuchebblocklanczos_eig(&ccm,&ccb);

  // print bands
  for(int ii=0; ii < nvecs; ii++){

    for(int jj=0; jj < nvecs+ccb.bsize; jj++){
      printf(" schurvecs[%d] = %+e\n", ii*(nvecs+ccb.bsize)+jj, ccb.schurvecs[ii*(nvecs+ccb.bsize)+jj]);
    }
    printf("\n");

    for(int jj=0; jj < ccb.bsize+1; jj++){
      printf(" bands[%d] = %+e\n", ii*(ccb.bsize+1)+jj,
             ccb.bands[ii*(ccb.bsize+1)+jj]);
    }
    printf("\n");

  }
  printf("\n");

  // print evals
  for(int ii=0; ii < nvecs; ii++){
    printf(" evals[%d] = %+e\n", ii, ccb.evals[ii]);
  }
  printf("\n");

  // print ritz vectors
  for(int jj=0; jj < nvecs; jj++){
  for(int ii=0; ii < ccb.n; ii++){
    cudaMemcpy(&val,&(ccb.dvecs)[jj*ccb.n + ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dvecs[%d] = %+e\n", jj*ccb.n+ii, val);
  }
  printf("\n");
  }
  printf("\n");

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCL
  cuchebblocklanczos_destroy(&ccb);

  // return 
  return 0;

}
