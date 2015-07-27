
#include "gpusollib.h"

#include "../SRC/spmv_csr_vector_kernel.cu"

void compare_mv(matrix_t *mat) {

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float time = 0;
  double flops = 0.0;

  csr_t *csr = mat->d_csr; // take csr format of matrix A
  int n = csr->n; // take size of matrix A
  int nnz = mat->nnz;

  double *y = (double*) cuda_malloc(n*sizeof(double));
  double *x = (double*) cuda_malloc(n*sizeof(double));
  cudaMemset(x, 1.0, n*sizeof(double));

  /*-------- The HALFWARP case */
  int hwb = BLOCKDIM / HALFWARP; // Determine number of half-warps per block
  int gDim = min(MAXTHREADS / BLOCKDIM, (n+hwb-1) / hwb);
  int bDim = BLOCKDIM;

  cudaEventRecord(start);
for (int i=0; i<100; i++) {
  spmv_csr_half_vector_kernel<<<gDim, bDim>>>(n, csr->ia, csr->ja, csr->a, x, y, 0);
}
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("Time to perform the MV - 16 threads per row:  %3.1f ms \n", time);

  time = 0;

  /*-------- The single thread case */
  cudaEventRecord(start);
for (int i=0; i<100; i++) {
  spmv_csr_scalar_kernel<<<ceil(n/BLOCKDIM), BLOCKDIM>>>(n, csr->ia, csr->ja, csr->a, x, y);
}
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("Time to perform the MV - 1 thread per row:  %3.1f ms \n", time);

  time = 0;

  /*-------- The WARP case */
  cudaEventRecord(start);
for (int i=0; i<100; i++) {
  spmv_csr_vector_kernel<<<ceil(n/BLOCKDIM)*WARP, BLOCKDIM>>>(n, csr->ia, csr->ja, csr->a, x, y);
}
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("Time to perform the MV - 32 threads per row:  %3.1f ms \n", time);

  time = 0;

/*-------- The 8 threads case */
  cudaEventRecord(start);
for (int i=0; i<100; i++) {
  spmv_csr_quarter_vector_kernel<<<ceil(n/BLOCKDIM)*8, BLOCKDIM>>>(n, csr->ia, csr->ja, csr->a, x, y);
}
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  printf("Time to perform the MV - 8 threads per row:  %3.1f ms \n", time);

  time = 0;





}
