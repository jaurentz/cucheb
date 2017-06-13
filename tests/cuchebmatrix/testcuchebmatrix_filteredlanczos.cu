#include <cucheb.h>

/* driver */
int main(){

  // read in matrix and allocate memory
  string mtxfile("../matrices/SiH4.mtx");
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // set interval [alpha,beta]
  double alpha, beta;
  //alpha = 9.0e6; beta = 1.0e8;  // easy
  //alpha = 4.0e6; beta = 5.0e6;  // less easy
  alpha = 2.5e6; beta = 3.0e6;  // even less easy

  // call filtered lanczos for an [alpha,beta]
  cucheblanczos ccl;
  cuchebstats ccstats;
  cuchebmatrix_filteredlanczos(alpha, beta, 1, &ccm, &ccl, &ccstats);

  // print matrix
  cuchebmatrix_print(&ccm);

  // print statistics
  cuchebstats_print(&ccstats);

  // print eigenvalues
  for (int ii=0; ii<ccl.nconv; ii++) {
    printf(" %+e, %e\n",ccl.evals[ccl.index[ii]],ccl.res[ccl.index[ii]]);
  } printf("\n");

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCB
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
