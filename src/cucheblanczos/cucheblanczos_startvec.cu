#include <cucheb.h>

/* routine to create starting vector for cucheblanczos */
int cucheblanczos_startvec(cucheblanczos* ccl){

  // initial starting vector to be normalized all ones vector
  double val;
  val = 1.0/sqrt((double)(ccl->n));
  for(int ii=0; ii < ccl->n; ii++){
    cudaMemcpy(&(ccl->dvecs)[ii],&val,sizeof(double),cudaMemcpyHostToDevice);
  }

  // return  
  return 0;

}
