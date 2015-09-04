#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("../matrices/H2O.mtx");
  string mtxfile("../matrices/G2_circuit.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  //string mtxfile("../matrices/Stranke94.mtx");

  // cuchebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);
  cuchebmatrix_specint(&ccm);
  cuchebmatrix_print(&ccm);

  // filter polynomial
  double tau;
  tau = 100.0*(ccm.m);
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);
  cuchebpoly_gaussianfilter(ccm.a,ccm.b,0,tau,&ccp);
  //cuchebpoly_pointfilter(ccm.a,ccm.b,0,100,&ccp);
  cuchebpoly_print(&ccp);

  // cuchebblocklanczos
  cuchebblocklanczos ccb;
  cuchebblocklanczos_init(3, 100, &ccm, &ccb);

  // print CCB
  cuchebblocklanczos_print(&ccb);

  // set starting vector
  cuchebblocklanczos_startvecs(&ccb);

  // do arnoldi run
  cuchebblocklanczos_filteredarnoldi(&ccm,&ccp,&ccb);

/*
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
*/

  // compute ritz values
  cuchebblocklanczos_eig(&ccm,&ccb);

/*
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
*/

  // compute rayleigh quotients and residuals
  cuchebblocklanczos_rayleigh(&ccm,&ccb);

  // print ritz vectors
  for(int jj=0; jj < (ccb.bsize)*(ccb.nblocks); jj++){
    printf(" evals[%d] = %+e, res[%d] = %+e\n", jj, ccb.evals[jj], 
           jj, ccb.res[jj]);
  }
  printf("\n");

  // destroy CCP
  cuchebpoly_destroy(&ccp);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCL
  cuchebblocklanczos_destroy(&ccb);

  // return 
  return 0;

}
