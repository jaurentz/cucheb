#include <cucheb.h>

/* routine to initialize cucheblanczos object */
int cucheblanczos_init(int bsize, int nblocks, cuchebmatrix* ccm, cucheblanczos* ccl){

  // set dimensions
  ccl->n = ccm->m;
  ccl->bsize = min(max(1,bsize),MAX_BLOCK_SIZE);
  ccl->nblocks = min(max(1,nblocks),MAX_NUM_BLOCKS);

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

