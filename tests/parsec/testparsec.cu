#include <cucheb.h>

/* driver */
int main(){

  // compute variables
  cuchebmatrix ccm;
  cucheblanczos ccl;

  // variables to parse file
  string mtxfile("../matrices/SiH4.mtx");
  double lbnd, ubnd;
  int bsize;

  // initialize matrix
  cuchebmatrix_init(mtxfile, &ccm);

  // set interval and block size
  lbnd = -0.645;
  ubnd = -0.0053;
  bsize = 3;

  // start timer
  clock_t tick;
  tick = clock();

  // call filtered lanczos for an interval
  cuchebmatrix_filteredlanczos(lbnd, ubnd, bsize, &ccm, &ccl);

  // computation time
  printf("\ncomputation time = %e\n",(clock()-tick)/((double)CLOCKS_PER_SEC));

  // print matrix
  cuchebmatrix_print(&ccm);

  // print eigenvalues and residuals
  printf("\nComputed eigenvalues and residuals:\n");
  for(int ii=0;ii<ccl.nconv;ii++){
    printf(" eig[%d] = %+e, res[%d] = %e\n",
           ii,ccl.evals[ccl.index[ii]],ii,ccl.res[ccl.index[ii]]);
  }

  // print first 10 entries of first eigenvector 
  printf("\nFirst 10 entries of first eigenvector:\n");
  for(int ii=0;ii<10;ii++){
    printf(" vec[%d] = %+e\n",ii,ccl.vecs[ccl.index[ii]*ccl.n+ii]);
  }

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // destroy ccm
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
