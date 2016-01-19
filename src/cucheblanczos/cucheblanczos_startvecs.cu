#include <cucheb.h>

/* routine to create starting vector for cucheblanczos */
int cucheblanczos_startvecs(cucheblanczos* ccl){

  // initial starting vector to be normalized all ones vector
  double scl;
  scl = 1.0/sqrt(1.0*ccl->n);
  for(int ii=0; ii < ccl->n; ii++){
    cudaMemcpy(&(ccl->dvecs)[ii],&scl,sizeof(double),cudaMemcpyHostToDevice);
  }

  // all other starting vectors
  double one = 1.0/sqrt(2.0), mone = -1.0/sqrt(2.0), zero = 0.0;
  for(int jj=1; jj < ccl->bsize; jj++){
    for(int ii=0; ii < ccl->n; ii++){
      if (ii == 2*(jj-1)) {
        cudaMemcpy(&(ccl->dvecs)[jj*(ccl->n) + ii],&one,sizeof(double),cudaMemcpyHostToDevice);
      }
      else if (ii == 2*(jj-1)+1) {
        cudaMemcpy(&(ccl->dvecs)[jj*(ccl->n) + ii],&mone,sizeof(double),cudaMemcpyHostToDevice);
      }
      else {
        cudaMemcpy(&(ccl->dvecs)[jj*(ccl->n) + ii],&zero,sizeof(double),cudaMemcpyHostToDevice);
      }
    }
  }

  // return  
  return 0;

}
