#include <gpusollib.h>


__global__
void spmv_csr_scalar_kernel( const int num_rows , const int * ptr, const int * indices, const double * data, const double * x, double * y) {

// An inefficient MV, one thread per row

int row = blockDim.x * blockIdx.x + threadIdx.x ;
if( row < num_rows ){

  float dot = 0;
  int row_start = ptr [ row ];
  int row_end = ptr [ row +1];
  for (int jj = row_start ; jj < row_end ; jj ++)
     dot += data [ jj ] * x[ indices [ jj ]];
  y[ row ] += dot ;
  }


}

__global__
void spmv_csr_vector_kernel(const int num_rows, const int * ptr, const int * indices, const double * data, const double * x, double * y) {

// Thirty-two threads per row

__shared__ float vals[BLOCKDIM];
int thread_id = blockDim.x * blockIdx.x + threadIdx.x; // global thread index
int warp_id = thread_id / 32; // global warp index
int lane = thread_id & (32 - 1); // thread index within the warp
int row = warp_id; // one warp per row

if (row < num_rows) {

  int row_start = ptr[row];
  int row_end = ptr[row+1];

  // compute running sum per thread
  vals[threadIdx.x] = 0;

  for(int jj = row_start + lane; jj < row_end; jj += 32)
     vals[threadIdx.x] += data[jj] * x[indices[jj]];

  // parallel reduction in shared memory
  if (lane < 16) vals[threadIdx.x] += vals[threadIdx.x + 16];
  if (lane < 8) vals[threadIdx.x] += vals[threadIdx.x + 8];
  if (lane < 4) vals[threadIdx.x] += vals[threadIdx.x + 4];
  if (lane < 2) vals[threadIdx.x] += vals[threadIdx.x + 2];
  if (lane < 1) vals[threadIdx.x] += vals[threadIdx.x + 1];

  // first thread writes the result
  if (lane == 0)
     y[row] += vals[threadIdx.x];
  }

}

__global__
void spmv_csr_quarter_vector_kernel(const int num_rows, const int * ptr, const int * indices, const double * data, const double * x, double * y) {

// Eight threads per row

__shared__ float vals[BLOCKDIM];
int thread_id = blockDim.x * blockIdx.x + threadIdx.x; // global thread index
int warp_id = thread_id / 8; // global warp index
int lane = thread_id & (8 - 1); // thread index within the warp
int row = warp_id; // one warp per row

if (row < num_rows) {

  int row_start = ptr[row];
  int row_end = ptr[row+1];

  // compute running sum per thread
  vals[threadIdx.x] = 0;

  for(int jj = row_start + lane; jj < row_end; jj += 8)
     vals[threadIdx.x] += data[jj] * x[indices[jj]];

  // parallel reduction in shared memory
  if (lane < 4) vals[threadIdx.x] += vals[threadIdx.x + 4];
  if (lane < 2) vals[threadIdx.x] += vals[threadIdx.x + 2];
  if (lane < 1) vals[threadIdx.x] += vals[threadIdx.x + 1];

  // first thread writes the result
  if (lane == 0)
     y[row] += vals[threadIdx.x];
  }

}

__global__  
void spmv_csr_half_vector_kernel(int n, int *d_ia, int *d_ja, double *d_a, double *d_x, double *d_y, int neg) { 
/*------------------------------------------------------------* 
 *               CSR spmv-vector kernel  
 *           Half-Warp (16 threads) per row 
 *------------------------------------------------------------*/ 
  // num of half-warps 
  int nhw = gridDim.x*BLOCKDIM/HALFWARP; 
  // half warp id 
  int hwid = (blockIdx.x*BLOCKDIM+threadIdx.x)/HALFWARP; 
  // thread lane in each half warp 
  int lane = threadIdx.x & (HALFWARP-1); 
  // half warp lane in each block 
  int hwlane = threadIdx.x/HALFWARP; 
  // shared memory for partial result 
  volatile __shared__ double r[BLOCKDIM+8]; 
  volatile __shared__ int startend[BLOCKDIM/HALFWARP][2]; 
 
  for (int row = hwid; row < n; row += nhw) { 
    // row start and end point 
    if (lane < 2) 
      startend[hwlane][lane] = d_ia[row+lane]; 
    int p = startend[hwlane][0]; 
    int q = startend[hwlane][1]; 
    double sum = 0.0; 
    for (int i=p+lane; i<q; i+=HALFWARP) { 
      sum += d_a[i-1] * d_x[d_ja[i-1]];
    } 
    // parallel reduction 
    r[threadIdx.x] = sum; 
    r[threadIdx.x] = sum = sum + r[threadIdx.x+8]; 
    r[threadIdx.x] = sum = sum + r[threadIdx.x+4]; 
    r[threadIdx.x] = sum = sum + r[threadIdx.x+2]; 
    r[threadIdx.x] = sum = sum + r[threadIdx.x+1]; 
    if (lane == 0) { 
      if (neg) 
        d_y[row] = -r[threadIdx.x]; 
      else 
        d_y[row] = r[threadIdx.x]; 
    } 
  } 
}

