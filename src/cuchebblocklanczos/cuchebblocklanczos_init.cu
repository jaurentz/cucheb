#include <cucheb.h>

/* routine to initialize cuchebblocklanczos object */
int cuchebblocklanczos_init(int bsize, int nblocks, cuchebmatrix* ccm, cuchebblocklanczos* ccb){

  // set dimensions
  ccb->n = ccm->m;
  ccb->bsize = min(max(1,bsize),MAX_BLOCK_SIZE);
  ccb->nblocks = min(max(1,nblocks),MAX_NUM_BLOCKS);

  // allocate host memory
  int nvecs;
  nvecs = (ccb->bsize)*(ccb->nblocks);
  ccb->index = new int[nvecs];
  if (ccb->index == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccb->bands = new double[nvecs*(ccb->bsize + 1)];
  if (ccb->bands == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccb->evals = new double[nvecs];
  if (ccb->evals == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  ccb->schurvecs = new double[nvecs*(nvecs + ccb->bsize)];
  if (ccb->schurvecs == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // initialize index
  for (int ii=0; ii < nvecs; ii++) {
    (ccb->index)[ii] = ii;
  }

  // allocate device memory
  if(cudaMalloc(&(ccb->dtemp),(nvecs + ccb->bsize)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccb->dvecs),(ccb->n)*(nvecs + ccb->bsize)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccb->dschurvecs),nvecs*(nvecs + ccb->bsize)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // return  
  return 0;

}

