#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("../matrices/Stranke94.mtx");

  // cuchebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);
  cuchebmatrix_specint(&ccm);
  cuchebmatrix_print(&ccm);

  // filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);
  cuchebpoly_pointfilter(ccm.a,ccm.b,0,10,&ccp);
  cuchebpoly_print(&ccp);

  // cucheblanczos
  cucheblanczos ccl;
  cucheblanczos_init(1, MAX_NUM_VECS, &ccm, &ccl);

  // print CCB
  cucheblanczos_print(&ccl);

  // set starting vector
  cucheblanczos_startvecs(&ccl);

  // cuchebstats
  cuchebstats ccstats;

  // do arnoldi run
  cucheblanczos_filteredarnoldi(5,&ccm,&ccp,&ccl,&ccstats);

  // print ccl
  cucheblanczos_print(&ccl);

  // print arnoldi vectors
  double val;
  int nvecs;
  nvecs = (ccl.bsize)*(ccl.nblocks);
  for(int jj=0; jj < nvecs+ccl.bsize; jj++){
  for(int ii=0; ii < ccl.n; ii++){
    cudaMemcpy(&val,&(ccl.dvecs)[jj*ccl.n + ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dvecs[%d] = %+e\n", jj*ccl.n+ii, val);
  }
  printf("\n");
  }
  printf("\n");

  // compute ritz values
  cucheblanczos_ritz(&ccm,&ccl);

  // print bands
  for(int ii=0; ii < nvecs; ii++){

    for(int jj=0; jj < nvecs+ccl.bsize; jj++){
      printf(" schurvecs[%d] = %+e\n", ii*(nvecs+ccl.bsize)+jj, ccl.schurvecs[ii*(nvecs+ccl.bsize)+jj]);
    }
    printf("\n");

    for(int jj=0; jj < ccl.bsize+1; jj++){
      printf(" bands[%d] = %+e\n", ii*(ccl.bsize+1)+jj,
             ccl.bands[ii*(ccl.bsize+1)+jj]);
    }
    printf("\n");

  }
  printf("\n");

  // print evals
  for(int ii=0; ii < nvecs; ii++){
    printf(" evals[%d] = %+e\n", ii, ccl.evals[ii]);
  }
  printf("\n");

  // print ritz vectors
  for(int jj=0; jj < nvecs; jj++){
  for(int ii=0; ii < ccl.n; ii++){
    cudaMemcpy(&val,&(ccl.dvecs)[jj*ccl.n + ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dvecs[%d] = %+e\n", jj*ccl.n+ii, val);
  }
  printf("\n");
  }
  printf("\n");

  // compute rayleigh quotients and residuals
  cucheblanczos_rayleigh(&ccm,&ccl);

  // print ritz vectors
  for(int jj=0; jj < (ccl.bsize)*(ccl.nblocks); jj++){
    printf(" evals[%d] = %+e, res[%d] = %+e\n", jj, ccl.evals[jj], 
           jj, ccl.res[jj]);
  }
  printf("\n");

  // destroy CCP
  cuchebpoly_destroy(&ccp);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCL
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
