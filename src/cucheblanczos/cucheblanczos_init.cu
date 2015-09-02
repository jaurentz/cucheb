#include <cucheb.h>

/* routine to initialize cucheblanczos object */
int cucheblanczos_init(int nvecs, cuchebmatrix* ccm, cucheblanczos* ccl){

  // set dimensions
  ccl->n = ccm->m;
  ccl->nvecs = min(ccl->n,max(nvecs,1));

  // allocate host memory
  ccl->diag = new double[ccl->nvecs];
  if (ccl->diag == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->sdiag = new double[ccl->nvecs];
  if (ccl->sdiag == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->schurvecs = new double[(ccl->nvecs)*(ccl->nvecs)];
  if (ccl->schurvecs == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // allocate device memory
  if(cudaMalloc(&(ccl->dtemp),(ccl->nvecs)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccl->dvecs),(ccl->n)*((ccl->nvecs)+1)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccl->dschurvecs),(ccl->nvecs)*(ccl->nvecs)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // return  
  return 0;

}

