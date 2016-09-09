#include <cucheb.h>
/*
  cucheblanczos_init

  This routine initializes an instance of a cucheblanczos object. The following
  inputs are required:

    bsize   - size of the Lanczos blocks
    numvecs - the total number of Lanczos vectors to allocate for
    ccm     - a reference to an initialized instance of a cuchebmatrix
    ccl     - a reference to an uninitialized instance of a cucheblanczos

*/

/* routine to initialize cucheblanczos object */
int cucheblanczos_init(int bsize, int numvecs, cuchebmatrix* ccm, cucheblanczos* ccl){

  // set dimensions
  ccl->n = ccm->m;
  ccl->bsize = min(max(1,bsize),MAX_BLOCK_SIZE);
  ccl->nblocks = min(min(ccl->n,MAX_NUM_VECS),max(1,numvecs))/(ccl->bsize);
  ccl->stop = 0;

  // allocate host memory
  int nvecs;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  ccl->index = new int[nvecs];
  if (ccl->index == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->evals = new double[nvecs];
  if (ccl->evals == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->res = new double[nvecs];
  if (ccl->res == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->bands = new double[nvecs*(ccl->bsize + 1)];
  if (ccl->bands == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->vecs = new double[nvecs*(ccl->n)];
  if (ccl->vecs == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccl->schurvecs = new double[nvecs*(nvecs + ccl->bsize)];
  if (ccl->schurvecs == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // initialize index
  for (int ii=0; ii < nvecs; ii++) {
    (ccl->index)[ii] = ii;
  }

  // allocate device memory
  if(cudaMalloc(&(ccl->dtemp),(nvecs + ccl->bsize)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccl->dvecs),(ccl->n)*(nvecs + ccl->bsize)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccl->dschurvecs),nvecs*(nvecs + ccl->bsize)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // return  
  return 0;

}

