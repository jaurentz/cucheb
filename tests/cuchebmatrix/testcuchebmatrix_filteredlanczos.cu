#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("../matrices/H2O.mtx");
  //string mtxfile("../matrices/Ga41As41H72.mtx");
  string mtxfile("../matrices/Si87H76.mtx");
  //string mtxfile("../matrices/Si34H36.mtx");
  //string mtxfile("../matrices/dielFilterV2real.mtx");
  //string mtxfile("../matrices/CO.mtx");
  //string mtxfile("../matrices/Si10H16.mtx");
  //string mtxfile("../matrices/G2_circuit.mtx");
  //string mtxfile("../matrices/Trefethen_20000.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // cucheblanczos
  cucheblanczos ccl;

time_t start = time(0);

  // call filtered lanczos for a point
  cuchebmatrix_filteredlanczos(4,-1e100,3,&ccm,&ccl);

  // call filtered lanczos for an interval
 // cuchebmatrix_filteredlanczos(-10.0, -1.0, 3, &ccm, &ccl);

time_t end = time(0);
double time = difftime(end, start);
printf(" \ntime = %e\n\n", time);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCB
  cucheblanczos_destroy(&ccl);

  // return 
  return 0;

}
