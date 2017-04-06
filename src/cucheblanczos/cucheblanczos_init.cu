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
    printf("Host memory allocation failed: index\n");
    exit(1);
  }
  ccl->evals = new double[nvecs];
  if (ccl->evals == NULL) {
    printf("Host memory allocation failed: evals\n");
    exit(1);
  }
  ccl->res = new double[nvecs];
  if (ccl->res == NULL) {
    printf("Host memory allocation failed: res\n");
    exit(1);
  }
  ccl->bands = new double[nvecs*(ccl->bsize + 1)];
  if (ccl->bands == NULL) {
    printf("Host memory allocation failed: bands\n");
    exit(1);
  }
  ccl->vecs = new double[nvecs*(ccl->n)];
  if (ccl->vecs == NULL) {
    printf("Host memory allocation failed: vecs\n");
    exit(1);
  }
  ccl->schurvecs = new double[nvecs*(nvecs + ccl->bsize)];
  if (ccl->schurvecs == NULL) {
    printf("Host memory allocation failed: schurvecs\n");
    exit(1);
  }

  // initialize index
  for (int ii=0; ii < nvecs; ii++) {
    (ccl->index)[ii] = ii;
  }

size_t freeMem, totalMem;
cudaMemGetInfo(&freeMem, &totalMem);
printf("cucheblanczos_init\n");
printf("Free = %ld, Total = %ld\n", freeMem, totalMem);

  // allocate device memory
  if(cudaMalloc(&(ccl->dtemp),(nvecs + ccl->bsize)*sizeof(double)) != 0) {
    printf("Device memory allocation failed: dtemp\n");
    exit(1);
  }
  if(cudaMalloc(&(ccl->dvecs),(ccl->n)*(nvecs + ccl->bsize)*sizeof(double)) != 0) {
    printf("Device memory allocation failed: dvecs\n");
    exit(1);
  }
  if(cudaMalloc(&(ccl->dschurvecs),nvecs*(nvecs + ccl->bsize)*sizeof(double)) != 0) {
    printf("Device memory allocation failed: dschurvecs\n");
    exit(1);
  }

cudaMemGetInfo(&freeMem, &totalMem);
printf("Free = %ld, Total = %ld\n", freeMem, totalMem);

  // return  
  return 0;

}

