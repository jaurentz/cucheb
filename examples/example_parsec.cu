#include <cucheb.h>

/* driver */
int main(){

  // compute variables
  cuchebmatrix ccm;
  cucheblanczos ccl;

  // variables to parse file
  string matname;
  double lbnd, ubnd;
  int bsize;

  // initialize matrix
  matname = "Ge87H76.mtx";
  cuchebmatrix_init(matname, &ccm);
  cuchebmatrix_print(&ccm);

  // set interval and block size
  lbnd = -0.645;
  ubnd = -0.0053;
  bsize = 3;

  // call filtered lanczos for an interval
  cuchebmatrix_filteredlanczos(lbnd, ubnd, bsize, &ccm, &ccl);

  // print eigenvalues and residuals
  printf("\nComputed eigenvalues and residuals:\n")
  for(int ii=0;ii<ccl.nconv;ii++){
    printf(" eig[%d] = %+e, res[%d] = %e\n",
           ii,ccl.evals[ccl.index[ii]],ii,ccl.res[ccl.index[ii]]);
  }

  // print first 10 entries of first eigenvector 
  printf("\nFirst 10 entries of first eigenvector:\n")
  for(int ii=0;ii<10;ii++){
    printf(" vec[%d] = %+e\n",ii,ccl.vecs[ccl.index[ii]*n+ii]);
  }

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // destroy ccm
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
