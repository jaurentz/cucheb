#include <gpusollib.h>

void filtered_spmv_csr_vector(matrix_t *mat, double *x, double *y, int degree, double w, double c, double *d_mu, int mat_choice) {

  // Initialization and declaration  
  csr_t  *csr = mat->d_csr;  // take csr format of matrix A
  int     n   = csr->n;      // take size of matrix A
  int     nnz = csr->nnz;    // take nonzero entries of matrix A
  double  scal, one = 1.0, zero = 0.0, minusc = 0.0;
  double  *vkm1, *vk, *Avk;
  double  minusone = -1.0;
  int     i, k;

  minusc = -c;

  cusparseMatDescr_t descra;
  cusparseCreateMatDescr(&descra);
  cusparseSetMatType(descra, CUSPARSE_MATRIX_TYPE_GENERAL);

  // For the CUSPARSE, CUBLAS contexts
  cusparseHandle_t cusparseHandle = 0;
  cusparseStatus_t cusparseStatus;
  cusparseStatus = cusparseCreate(&cusparseHandle);

  /*
  if (checkCudaErrors(cusparseStatus))
  {
     exit(EXIT_FAILURE);
  }
  */

  cublasHandle_t cublasHandle = 0;
  cublasStatus_t cublasStatus;
  cublasStatus = cublasCreate(&cublasHandle);

  /*
  if (checkCudaErrors(cublasStatus))
  {
    exit(EXIT_FAILURE);
  }
  */

  /*
  int hwb  = BLOCKDIM / HALFWARP; // Determine number of half-warps per block
  int gDim = min(MAXTHREADS / BLOCKDIM, (n+hwb-1) / hwb);
  int bDim = BLOCKDIM;
  */


  // allocate memory
  cudaMalloc( (double**)&vkm1, n*sizeof(double) );
  cudaMalloc( (double**)&vk,   n*sizeof(double) );
  cudaMalloc( (double**)&Avk,  n*sizeof(double) );

  // initialize buffers
  cudaMemset( vkm1, 0.0, n*sizeof(double) );
  cudaMemset( vk, 0.0,   n*sizeof(double) );
  cudaMemset( Avk, 0.0,  n*sizeof(double) );
  cudaMemset( y, 0.0,    n*sizeof(double) );

  // copy x to vk
  cudaMemcpy( vk, x, n*sizeof(double), cudaMemcpyDeviceToDevice );


double *temp = (double*) malloc((nnz)*sizeof(double));


  for (k = 0; k <= degree; k++) {
   
     cublasDaxpy(cublasHandle, n, &d_mu[k], vk, 1, y, 1);

     scal = 2.0 / w;
     if (k==0)
        scal = 1.0 / w;

     if ( mat_choice == 0 ) {
        cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, nnz, &one, descra, csr->a, csr->ia, csr->ja, vk, &zero, Avk);
     }

     cublasDaxpy(cublasHandle, n, &minusc, vk, 1, Avk, 1);
     cublasDscal(cublasHandle, n, &scal, Avk, 1);
     cublasDaxpy(cublasHandle, n, &minusone, vkm1, 1, Avk, 1);

     cudaMemcpy( vkm1, vk, n*sizeof(double), cudaMemcpyDeviceToDevice );  // vkm1 = vk;
     cudaMemcpy( vk, Avk, n*sizeof(double), cudaMemcpyDeviceToDevice  );  // vk = vkp1;
     cudaMemset( Avk, 0.0, n*sizeof(double) );

  }

  cusparseDestroy(cusparseHandle);
  cublasDestroy(cublasHandle);
  cudaFree(vk);
  cudaFree(vkm1);
  cudaFree(Avk);

}










