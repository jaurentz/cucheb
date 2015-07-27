#include <gpusollib.h>
#include <texture.h>

__global__ 
void csr_v_k(int n, int *d_ia, int *d_ja, double *d_a, 
             double *d_y, int neg) {
/*------------------------------------------------------------*
 *               CSR spmv-vector kernel
 *  shared memory reduction, texture memory fetching
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
      int2 t = tex1Dfetch(tex_double, d_ja[i-1]-1);
      sum += d_a[i-1] * __hiloint2double(t.y, t.x);
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

/*-----------------------------------------------*/
void spmv_csr_vector(matrix_t *mat, double *x, 
                     double *y, int neg) {
/*-----------------------------------------------*/
  csr_t *csr = mat->d_csr;
  int n = csr->n;
/*--------- texture binding */
  size_t offset;
  cuda_bind_tex(&offset, x, n*sizeof(int2));
  assert(offset == 0);
/*-------- set spmv kernel */
/*-------- num of half-warps per block */
  int hwb = BLOCKDIM/HALFWARP;
  int gDim = min(MAXTHREADS/BLOCKDIM, (n+hwb-1)/hwb);
  int bDim = BLOCKDIM;
/*-------------*/
  csr_v_k<<<gDim, bDim>>>
  (n, csr->ia, csr->ja, csr->a, y, neg);
/*--------- unbind texture */
  cuda_unbind_tex();
}

/*----------------------------------------------*/
/*             THE JAD FORMAT                   */
/*----------------------------------------------*/
__global__ 
void jad_k(int n, int njad, int *d_ia, int *d_ja, 
           double *d_a, double *d_y, int neg) {
/*-----------------------------------------------*
 *  JAD SpMatVec Kernel(texture memory fetching)
 *            one thread per row
 *-----------------------------------------------*/
  int i,j,p,q;
  double r;
/*------------ each thread per row */
  int row = blockIdx.x*blockDim.x+threadIdx.x;
/*------------ number of threads */
  int nthreads = gridDim.x * blockDim.x;
  __shared__ int shia[BLOCKDIM];
  if (threadIdx.x <= njad)
    shia[threadIdx.x] = d_ia[threadIdx.x];
  __syncthreads(); 
  while (row < n) {
    r = 0.0;
    p = shia[0];
    q = shia[1];
    i=0;
    while ( ((p+row) < q) && (i < njad) ) {
      j = p+row;
/*--------- double precision texture fetching */
      int2 t = tex1Dfetch(tex_double, d_ja[j-1]-1);
      r += d_a[j-1] * __hiloint2double(t.y, t.x);
      i++;
      if (i<njad) {
        p = q;
        q = shia[i+1];
      }
    }
    if (neg)
      d_y[row] = -r;
    else
      d_y[row] = r;
    row += nthreads;
  }
}

__global__ 
void vperm_k(int n, double *d_y, double *d_x, int *d_p) {
/*------------------------------------------------*/
/*   vector permutation, y[p(i)] := x[i]          */
/*------------------------------------------------*/
  int idx = blockIdx.x*blockDim.x+threadIdx.x;
  int nthreads = gridDim.x * blockDim.x;
  int i;
  for (i=idx; i<n; i+=nthreads)
    d_y[d_p[i]-1] = d_x[i];
}

/*------------------------------------------*/
void spmv_jad(matrix_t *mat, double *x, 
              double *y, int neg) {
/*------------------------------------------*/
  jad_t *jad = mat->d_jad;
  int n = jad->n;
/*--------- texture binding */
  size_t offset;
  cuda_bind_tex(&offset, x, n*sizeof(int2));
  assert(offset == 0);
/*--------- set spmv kernel */
  int nthreads = min(MAXTHREADS, n);
  int gDim = (nthreads+BLOCKDIM-1)/BLOCKDIM;
  int bDim = BLOCKDIM;
/*--------- w = A*x */
  jad_k<<<gDim, bDim>>>
  (n, jad->njad, jad->ia, jad->ja, jad->a, 
   jad->w, neg);
/*-------- permutation */
  vperm_k<<<gDim, bDim>>>
  (n, y, jad->w, jad->perm);
/*--------- unbind texture */
  cuda_unbind_tex();
}

/*---------------------------------------------------*/
/*                THE DIA FORMAT                     */
/*---------------------------------------------------*/
__global__ 
void dia_k(int n, int stride, int ndiags, int *d_ioff, 
              double *d_diags, double *d_y, int neg) {
/*---------------------------------------------------*
 *    DIA SpMatVec Kernel(texture memory fetching)
 *               one thread per row
 *---------------------------------------------------*/
/*--------- number of threads */
  int nthreads = gridDim.x * blockDim.x;
  __shared__ int shioff[BLOCKDIM];
  if (threadIdx.x < ndiags)
    shioff[threadIdx.x] = d_ioff[threadIdx.x];
  __syncthreads();
/*--------- each thread per row */
  int row = blockIdx.x*blockDim.x+threadIdx.x;
/*--------- start position of a row */
  d_diags += row;
  while (row < n) {
    double r = 0;
    for (int i=0; i<ndiags; i++) {
/*--------- col index */
      int q = shioff[i] + row;
      if (q >= 0 && q < n) {
/*--------- double precision texture fetching */
        int2 t = tex1Dfetch(tex_double, q);
        r += *(d_diags + i*stride) * __hiloint2double(t.y, t.x);
      }
    }
    if (neg)
      d_y[row] = -r;
    else
      d_y[row] = r;
    row += nthreads;
    d_diags += nthreads;
  }
}

/*---------------------------------------*/
void spmv_dia(matrix_t *mat, double *x, 
double *y, int neg) {
/*---------------------------------------*/
  dia_t *dia = mat->d_dia;
/*-------- GPU DIA format SpMV */
  int n = dia->n;
/*--------- texture binding */
  size_t offset;
  cuda_bind_tex(&offset, x, n*sizeof(int2));
  assert(offset == 0);
/*------- set spmv kernel */
  int nthreads = min(MAXTHREADS, n);
  int gDim = (nthreads+BLOCKDIM-1)/BLOCKDIM;
  int bDim = BLOCKDIM;
/*-------- */
  dia_k<<<gDim, bDim>>>
  (n, dia->stride, dia->ndiags, dia->ioff,
   dia->diags, y, neg);
/*--------- unbind texture */
  cuda_unbind_tex();
}

/*-------------------------------------------------*/
void spmv_csr_cpu(csr_t *csr, double *x, double *y) {
/*------------- CPU CSR SpMV kernel */
  int i,j;
  double r;
  for (i=0; i<csr->n; i++) {
    r = 0.0;
    for (j=csr->ia[i]; j<csr->ia[i+1]; j++)
      r += csr->a[j-1]*x[csr->ja[j-1]-1];
    y[i] = r;
  }
}

/*--------------------*/
/*  SpMV for MC-SSOR  */
/*  w = A*y-x         */
/*--------------------*/
__global__ 
void spmv_sor_k(int n, int *d_ia, int *d_ja, 
double *d_a,  double *d_x, double *d_w) {
/*----------------------------------------*/
  // num of half-warps
  int nhw = gridDim.x*BLOCKDIM/HALFWARP;
  // half warp id
  int hwid = (blockIdx.x*BLOCKDIM+threadIdx.x)/HALFWARP;
  // thread lane in each half warp
  int lane = threadIdx.x & (HALFWARP-1);
  // half warp lane in each block
  int hwlane = threadIdx.x/HALFWARP;
  // shared memory for patial result
  volatile __shared__ double r[BLOCKDIM+8];
  volatile __shared__ int startend[BLOCKDIM/HALFWARP][2];
/*------------------------------------------*/
  for (int row = hwid; row < n; row += nhw) {
    // row start and end point
    if (lane < 2)
      startend[hwlane][lane] = d_ia[row+lane];
    int p = startend[hwlane][0];
    int q = startend[hwlane][1];
    double sum = 0.0;
    for (int i=p+lane; i<q; i+=HALFWARP) {
      int2 t = tex1Dfetch(tex_double, d_ja[i-1]-1);
      sum += d_a[i-1] * __hiloint2double(t.y, t.x);
    }
    // parallel reduction
    r[threadIdx.x] = sum;
    r[threadIdx.x] = sum = sum + r[threadIdx.x+8];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+4];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+2];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+1];
    if (lane == 0)
      d_w[row] = d_x[row] - r[threadIdx.x];
  }
}

/*----------------------------------------------------*/
void spmv_sor(int n, int nrow, int *d_ia, int *d_ja, 
double *d_a, double *d_y, double *d_x, double *d_w) {
/*----------------------------------------------------*/
/*--------- texture binding */
  size_t offset;
  cuda_bind_tex(&offset, d_y, n*sizeof(int2));
  assert(offset == 0);
/*------ set Spmv kernel */
/*------ num of half-warps per block */
  int hwb = BLOCKDIM/HALFWARP; 
  int gDim = min(MAXTHREADS/BLOCKDIM, (nrow+hwb-1)/hwb);
  int bDim = BLOCKDIM;
/*---------------*/
  spmv_sor_k<<<gDim, bDim>>>
  (nrow, d_ia, d_ja, d_a, d_x, d_w);
/*--------- unbind texture */
  cuda_unbind_tex();
}

/*---------------*/
/* Spmv for ILU0 */
/*---------------*/
__global__ 
void spmv_ilu0_k1(int n, int *d_ia, int *d_ja, 
double *d_a, double *d_y) {
/*------------------------------------------*/
  // num of half-warps
  int nhw = gridDim.x*BLOCKDIM/HALFWARP;
  // half warp id
  int hwid = (blockIdx.x*BLOCKDIM+threadIdx.x)/HALFWARP;
  // thread lane in each half warp
  int lane = threadIdx.x & (HALFWARP-1);
  // half warp lane in each block
  int hwlane = threadIdx.x/HALFWARP;
  // shared memory for patial result
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
      int2 t = tex1Dfetch(tex_double, d_ja[i-1]-1);
      sum += d_a[i-1] * __hiloint2double(t.y, t.x);
    }
    // parallel reduction
    r[threadIdx.x] = sum;
    r[threadIdx.x] = sum = sum + r[threadIdx.x+8];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+4];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+2];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+1];
    if (lane == 0)
      d_y[row] -= r[threadIdx.x];
  }
}

/*-----------------------------------------------*/
/* d_y(noff:noff+nrow) -= A(noff:noff+nrow) * d_y*/
/*-----------------------------------------------*/
void spmv_ilu0_1(int noff, int nrow, int *d_ia, 
int *d_ja, double *d_a, double *d_y) {
/*-----------------------------------------------*/
  size_t offset;
  cuda_bind_tex(&offset, d_y, noff*sizeof(int2));
  assert(offset == 0); 
/*------- set spmv kernel */
  int hwb = BLOCKDIM/HALFWARP;
  int gDim = min(MAXTHREADS/BLOCKDIM, (nrow+hwb-1)/hwb);
  int bDim = BLOCKDIM;
/*------------------*/
  spmv_ilu0_k1<<<gDim, bDim>>>
  (nrow, d_ia, d_ja, d_a, d_y+noff);
/*------------------*/
  cuda_unbind_tex();
}

/*-------------------------------------------*/
__global__ 
void spmv_ilu0_k2(int n, int *d_ia, int *d_ja, 
double *d_a, double *d_w) {
/*-------------------------------------------*/
  // num of half-warps
  int nhw = gridDim.x*BLOCKDIM/HALFWARP;
  // half warp id
  int hwid = (blockIdx.x*BLOCKDIM+threadIdx.x)/HALFWARP;
  // thread lane in each half warp
  int lane = threadIdx.x & (HALFWARP-1);
  // half warp lane in each block
  int hwlane = threadIdx.x/HALFWARP;
  // shared memory for patial result
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
      int2 t = tex1Dfetch(tex_double, d_ja[i-1]-1);
      sum += d_a[i-1] * __hiloint2double(t.y, t.x);
    }
    // parallel reduction
    r[threadIdx.x] = sum;
    r[threadIdx.x] = sum = sum + r[threadIdx.x+8];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+4];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+2];
    r[threadIdx.x] = sum = sum + r[threadIdx.x+1];
    if (lane == 0)
      d_w[row] = -r[threadIdx.x];
  }
}

/*--------------------------------------------------*/
void spmv_ilu0_2(int n, int nrow, int *d_ia, int *d_ja, 
double *d_a, double *d_y, double *d_w) {
/*--------------------------------------------------*/
  size_t offset;
  cuda_bind_tex(&offset, d_y, n*sizeof(int2));
  assert(offset == 0); 
/*------- set spmv kernel */
  int hwb = BLOCKDIM/HALFWARP;
  int gDim = min(MAXTHREADS/BLOCKDIM, (nrow+hwb-1)/hwb);
  int bDim = BLOCKDIM;
  spmv_ilu0_k2<<<gDim, bDim>>>
  (nrow, d_ia, d_ja, d_a, d_w);
/*---------------------- */
  cuda_unbind_tex();
}

