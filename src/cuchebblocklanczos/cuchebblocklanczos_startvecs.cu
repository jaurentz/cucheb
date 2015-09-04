#include <cucheb.h>

/* routine to create starting vector for cuchebblocklanczos */
int cuchebblocklanczos_startvecs(cuchebblocklanczos* ccb){

  // initial starting vector to be normalized all ones vector
  double one = 1.0, zero = 0.0;
  for(int jj=0; jj < ccb->bsize; jj++){
    for(int ii=0; ii < ccb->n; ii++){
      if (ii == jj) {
        cudaMemcpy(&(ccb->dvecs)[jj*(ccb->n) + ii],&one,sizeof(double),cudaMemcpyHostToDevice);
      }
      else {
        cudaMemcpy(&(ccb->dvecs)[jj*(ccb->n) + ii],&zero,sizeof(double),cudaMemcpyHostToDevice);
      }
    }
  }

  // return  
  return 0;

}
