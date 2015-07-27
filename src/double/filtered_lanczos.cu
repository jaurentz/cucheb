
#include "gpusollib.h"
#define  SEED 100


void filtered_lanczos(matrix_t *mat, int msteps, int degree, double w, double c, int mat_choice) {


   //-------------------
   // Declare variables
   //------------------
   int     i, n;
   double  *h_v, *d_w, *d_u, *VV, *d, *e, *z, *work, *temp;
   double  t, dummy, beta, minusbeta, alpha, minusalpha, orthTol, wn, one = 1.0, zero = 0.0, minusone = -1.0;
   double  alpha2, beta2, t2, dummy2, minusalpha2, minusbeta2, d_one, d_zero, d_minusone;

   cublasHandle_t cublasHandle = 0;
   cublasStatus_t cublasStatus;
   cublasStatus = cublasCreate(&cublasHandle);

   //-----------------------------------
   // Initialize random number generator
   //-----------------------------------
   srand(SEED);

   // Number of rows
   n = mat->n;

   // Generate initial random vector for Lanczos
   h_v = (double*) malloc( n * sizeof(double) );
   for (i = 0; i < n; i++) {
      h_v[i] = rand()/((double)RAND_MAX + 1);
   }

   //---------------------------------------------------------
   // Set up device memory for Lanczos basis and other vectors
   //---------------------------------------------------------
   cudaMalloc( (void**)&VV, (msteps+1) * n * sizeof(double) );
   cudaMemset(VV, 0.0, (msteps+1) * n * sizeof(double) );

   cudaMalloc( (void**)&d_w, n * sizeof(double) );
   cudaMemset(d_w, 0.0, n * sizeof(double) );

   cudaMalloc( (void**)&d_u, msteps * sizeof(double) );
   cudaMemset(d_u, 0.0, n * sizeof(double) );

   cudaMalloc( (void**)&temp, n * sizeof(double) );
   cudaMemset(temp, 0.0, n * sizeof(double) );

   // Scalar values that reside in the GPU
   cudaMalloc( (void**)&alpha2, sizeof(double) );
   cudaMemset(&alpha2, 0.0, sizeof(double) );
   cudaMalloc( (void**)&beta2, sizeof(double) );
   cudaMemset(&beta2, 0.0, sizeof(double) );
   cudaMalloc( (void**)&t2, sizeof(double) );
   cudaMemset(&t2, 0.0, sizeof(double) );
   cudaMalloc( (void**)&minusalpha2, sizeof(double) );
   cudaMemset(&minusalpha2, 0.0, sizeof(double) );
   cudaMalloc( (void**)&minusbeta2, sizeof(double) );
   cudaMemset(&minusbeta2, 0.0, sizeof(double) );
   cudaMalloc( (void**)&dummy2, sizeof(double) );
   cudaMemset(&dummy2, 0.0, sizeof(double) );
   cudaMalloc( (void**)&d_one, sizeof(double) );
   cudaMemset(&d_one, 1.0, sizeof(double) );
   cudaMalloc( (void**)&d_zero, sizeof(double) );
   cudaMemset(&d_zero, 0.0, sizeof(double) );
   cudaMalloc( (void**)&d_minusone, sizeof(double) );
   cudaMemset(&d_minusone, -1.0, sizeof(double) );

   // Scale starting Lanczos vector
   cudaMemcpy( VV, h_v, n * sizeof(double), cudaMemcpyHostToDevice );
   cublasDnrm2(cublasHandle, n, VV, 1, &t2);
   cudaMemcpy( &t, &t2, n * sizeof(double), cudaMemcpyDeviceToHost );
   dummy = 1.0 / t;
   cudaMemcpy( &dummy2, &dummy, n * sizeof(double), cudaMemcpyHostToDevice );
   cublasDscal(cublasHandle, n, &dummy2, VV, 1);

   // set up other initial variables
   beta    = 0.0;
   alpha   = 0.0;
   orthTol = 1.0e-8; // just some checking
   wn      = 0.0;    //    >>


   //------------------------------------------------------
   // Allocate space to hold the tridiagonal Lanczos matrix
   // -----------------------------------------------------
   d = (double*) malloc( msteps * sizeof(double) );
   e = (double*) malloc( (msteps-1) * sizeof(double) );


   //---------------------------------
   // Description of filter parameters
   //---------------------------------
   // damping  : Jackson, Jackson-Chebychev, Delta 
   // xi       : test vector to check the quality of the polynomial
   // [i1, i2] : After scaling the eigvals of A, we seek the 
   //            eigenvalues in [i1,i2]
   // mu       : the coefficients
   // d_mu     : the coefficients' copy in the device

   int     damping;
   double  i1, i2;
   double  *mu, *d_mu;
   i1 = 0.15;
   i2 = 0.20;
   damping = 0;
   mu = (double*) malloc( (degree + 1) * sizeof(double) );
   memset (mu, 0.0, (degree + 1) * sizeof(double) );


   //--------------------------------------------------------
   // Compute polynomial coefficients and copy them to device
   compute_coeff(degree, i1, i2, damping, mu);
   cudaMalloc( (void**)&d_mu, (degree+1) * sizeof(double) );
   cudaMemcpy( d_mu, mu, (degree+1)*sizeof(double), cudaMemcpyHostToDevice );
  
   //--------------------------
   // Lanczos phase has started
   //--------------------------
   printf("Lanczos Alg begins ...\n");

   for ( i = 0; i < msteps; i++) {
     //-------------------
     // MV -- d_w = A*temp
     //-------------------
     cudaMemcpy(temp, &VV[i*n], n * sizeof(double), cudaMemcpyDeviceToDevice);
     filtered_spmv_csr_vector(mat, temp, d_w, 0, degree, w, c, d_mu, mat_choice); 

     //------------------
     // 3-term recurrence
     //------------------
     minusbeta = -beta;
     cudaMemcpy( &minusbeta2, &minusbeta, sizeof(double), cudaMemcpyHostToDevice );
     if (i == 0)
        cublasDaxpy(cublasHandle, n, &minusbeta2, VV, 1, d_w, 1);
     else
        cublasDaxpy(cublasHandle, n, &minusbeta2, &VV[(i-1)*n], 1, d_w, 1);

     // Compute alpha
     cublasDdot (cublasHandle, n, d_w, 1, &VV[i*n], 1, &alpha2);
     cudaMemcpy( &alpha, &alpha2, sizeof(double), cudaMemcpyDeviceToHost );
     minusalpha = -alpha;
     cudaMemcpy( &minusalpha2, &minusalpha, sizeof(double), cudaMemcpyHostToDevice );
     cublasDaxpy(cublasHandle, n, &minusalpha2, &VV[i*n], 1, d_w, 1);

     // Add on-diagonal entry
     d[i] = alpha;
     wn += alpha*alpha;
 
     //---------------------------
     // Re-orthogonalization phase
     //---------------------------
     // u = V'*w
     cublasDgemv(cublasHandle, CUBLAS_OP_T, n, i+1, &d_one, VV, n, d_w, 1, &d_zero, d_u, 1);

     // w = w - V*u
     cublasDgemv(cublasHandle, CUBLAS_OP_N, n, i+1, &d_minusone, VV, n, d_u, 1, &d_one, d_w, 1);

     //-----------------------------------------------
     // Take norm of vector after re-orthogonalization
     //-----------------------------------------------
     cublasDdot (cublasHandle, n, d_w, 1, d_w, 1, &beta2);
     cudaMemcpy( &beta, &beta2, sizeof(double), cudaMemcpyDeviceToHost );

     if (beta*(i+1) < orthTol*wn)
       break;
     wn += 2.0 * beta;

     //-----------------
     // Normalize vector
     //-----------------
     beta = sqrt(beta);
     dummy = 1.0 / beta;
     cudaMemcpy( &dummy2, &dummy, sizeof(double), cudaMemcpyHostToDevice );
     cublasDaxpy(cublasHandle, n, &dummy2, d_w, 1, &VV[(i+1)*n], 1);

     // Add new off-diagonal entry to T
     if (i < msteps-1) {
       e[i] = beta;
     }

   }


   //----------------------
   // Lanczos phase is over
   //----------------------
   printf("Generating the Lanczos basis : Done\n");


   //-----------------------------------------
   // Solve the tridiagonal eigenvalue problem
   //-----------------------------------------
   printf("Lapack STEQR begins ...\n");
   work = (double*) malloc( 2 * (msteps-1) * sizeof(double) );
   Calloc(z, msteps*msteps, double);
   int info;
   char compz = 'I';
   STEQR(&compz, &msteps, d, e, z, &msteps, work, &info); // use BLAS routines
   if (info != 0) {
     printf("LAPACK: FAILED TO FIND EIGENVALUES !!!\n");
     exit(-1);
   }


   //-------------------
   // Deallocate vectors
   //-------------------
   cudaFree(VV);
   cudaFree(d_mu);
   cudaFree(d_w);
   cudaFree(d_u);
   cudaFree(temp);
   free(h_v);
   free(d);
   free(e);
   free(work);
   free(z);
   free(mu);


}






