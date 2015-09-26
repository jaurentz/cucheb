#include <cucheb.h>

/* driver */
int main(){

  // set device
  //cudaSetDevice(1);

  // input file
  //string mtxfile("../matrices/SiH4.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  //string mtxfile("../matrices/H2O.mtx");
  //string mtxfile("../matrices/Si34H36.mtx");
  //string mtxfile("../matrices/Si87H76.mtx");
  //string mtxfile("../matrices/CO.mtx");
  //string mtxfile("../matrices/Ga41As41H72.mtx");
  //string mtxfile("../matrices/Ge99H100.mtx");
  //string mtxfile("../matrices/Andrews.mtx");
  //string mtxfile("../matrices/Laplacian.mtx");
  //string mtxfile("../matrices/Qdot3.mtx");
  //string mtxfile("../matrices/DIMACS/144.mtx");
  //string mtxfile("../matrices/DIMACS/598a.mtx");
  string mtxfile("../matrices/DIMACS/al2010.mtx");
  //string mtxfile("../matrices/DIMACS/ar2010.mtx");
  //string mtxfile("../matrices/DIMACS/auto.mtx");
  //string mtxfile("../matrices/DIMACS/az2010.mtx");
  //string mtxfile("../matrices/DIMACS/ca2010.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cucheblanczos
  cucheblanczos ccl;

  // cuchebstats
  cuchebstats ccstats;

  //cuchebmatrix_specint(&ccm);

  // call filtered lanczos for a point
  //cuchebmatrix_filteredlanczos(100, 1.0e300, 1, &ccm, &ccl, &ccstats);

  // call filtered lanczos for an interval
  //cuchebmatrix_filteredlanczos(-0.66, -.33, 3, &ccm, &ccl, &ccstats);

  // call expert lanczos
  //cuchebmatrix_expertlanczos(4.00, 5.00, 150, 1, 4000, 4000, &ccm, &ccl, &ccstats);
  //cuchebmatrix_expertlanczos(1.00, 1.01, 1600, 3, 1200, 400, &ccm, &ccl, &ccstats);
  //cuchebmatrix_expertlanczos(-2.0, -.33, 100, 3, 1200, 40, &ccm, &ccl, &ccstats);
  cuchebmatrix_expertlanczos(1.0e6, 1.0e100, -1, 1, 3000, 500, &ccm, &ccl, &ccstats);

  // print ccm
  cuchebmatrix_print(&ccm);

  // print ccstats
  cuchebstats_print(&ccstats);

  // print eigenvalues
  for (int ii=0; ii<ccl.nconv; ii++) {
  //for (int ii=0; ii<100; ii++) {
    printf(" %+e, %e\n",ccl.evals[ccl.index[ii]],ccl.res[ccl.index[ii]]);
  }
  printf("\n");

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCB
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
