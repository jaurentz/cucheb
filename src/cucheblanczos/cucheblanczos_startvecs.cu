#include <cucheb.h>

/* routine to create starting vector for cucheblanczos */
int cucheblanczos_startvecs(cucheblanczos* ccl){

  // initial starting vector to be normalized all ones vector
  double one = 1.0, zero = 0.0;
  for(int jj=0; jj < ccl->bsize; jj++){
    for(int ii=0; ii < ccl->n; ii++){
      if (ii == jj) {
        cudaMemcpy(&(ccl->dvecs)[jj*(ccl->n) + ii],&one,sizeof(double),cudaMemcpyHostToDevice);
      }
      else {
        cudaMemcpy(&(ccl->dvecs)[jj*(ccl->n) + ii],&zero,sizeof(double),cudaMemcpyHostToDevice);
      }
    }
  }

  // return  
  return 0;

}
