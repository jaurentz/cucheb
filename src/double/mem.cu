#include "gpusollib.h"

/*-----------------------------------*/
void *cuda_malloc(int size) {
  void *p = NULL;
  cudaMalloc(&p, size);
  return(p);
}

/*---------------------------------------*/
void *cuda_malloc_host(int size) {
  void *p = NULL;
  cudaMallocHost(&p, size);
  return (p);
}

/*---------------------------------------------*/
void memcpyh2d(void *dest, void* src, int size) {
  cudaMemcpy(dest, src, size, cudaMemcpyHostToDevice);
}

/*---------------------------------------------*/
void memcpyd2h(void *dest, void *src, int size) {
  cudaMemcpy(dest, src, size, cudaMemcpyDeviceToHost);
}

/*---------------------------------------------*/
void memcpyd2d(void *dest, void *src, int size) {
  cudaMemcpy(dest, src, size,
  cudaMemcpyDeviceToDevice);
}

/*---------------------------------------------*/
void cuda_memset(void *addr, int val, int size) {
  cudaMemset(addr, val, size);
}

/*------------------------------------------*/
void malloc_csr(int n, int nnz, csr_t *csr) {
  csr->n = n;
  csr->nnz = nnz;
  Malloc(csr->ia, n+1, int);
  Malloc(csr->ja, nnz, int);
  Malloc(csr->a,  nnz, double);
}

/*-----------------------------------------*/
void realloc_csr(csr_t *csr, int nnz) {
  csr->nnz = nnz;
  Realloc(csr->ja, nnz, int);
  Realloc(csr->a, nnz, double);
}

/*-----------------------------------------*/
void cuda_malloc_csr(int n, int nnz, 
                     csr_t *d_csr) {
  d_csr->ia = (int*)cuda_malloc((n+1)*sizeof(int));
  d_csr->ja = (int*)cuda_malloc(nnz*sizeof(int));
  d_csr->a = (double*)cuda_malloc(nnz*sizeof(double));
}

/*-------------------------------------------*/
void copy_csr_h2d(csr_t *h_csr, csr_t *d_csr) {
  int n = h_csr->n;
  int nnz = h_csr->nnz;
  d_csr->n = n;
  d_csr->nnz = nnz;
  memcpyh2d(d_csr->ia, h_csr->ia, (n+1)*sizeof(int));
  memcpyh2d(d_csr->ja, h_csr->ja, nnz*sizeof(int));
  memcpyh2d(d_csr->a,  h_csr->a,  nnz*sizeof(double));
}

/*--------------------------------------------*/
void copy_csr_h2h(csr_t *csr1, csr_t *csr2) {
  int n = csr1->n;
  int nnz = csr1->nnz;
  csr2->n = n;
  csr2->nnz = nnz;
  memcpy(csr2->ia, csr1->ia, (n+1)*sizeof(int));
  memcpy(csr2->ja, csr1->ja, nnz*sizeof(int));
  memcpy(csr2->a,  csr1->a,  nnz*sizeof(double));
}

/*------------------------------------------*/
void cuda_malloc_jad(int n, int njad, 
                     int nnz, jad_t *d_jad) {
  d_jad->ia = (int*)cuda_malloc((njad+1)*sizeof(int));
  d_jad->ja = (int*)cuda_malloc(nnz*sizeof(int));
  d_jad->a  = (double*)cuda_malloc(nnz*sizeof(double));
  d_jad->perm = (int*)cuda_malloc(n*sizeof(int));
  d_jad->w = (double*)cuda_malloc(n*sizeof(double));
}

/*--------------------------------------------*/
void copy_jad_h2d(jad_t *h_jad, jad_t *d_jad) {
  int njad = h_jad->njad;
  int nnz = h_jad->nnz;
  int n = h_jad->n;
  d_jad->n = n;
  d_jad->nnz = nnz;
  d_jad->njad = njad;
  memcpyh2d(d_jad->ia, h_jad->ia, (njad+1)*sizeof(int));
  memcpyh2d(d_jad->ja, h_jad->ja, nnz*sizeof(int));
  memcpyh2d(d_jad->a, h_jad->a, nnz*sizeof(double));
  memcpyh2d(d_jad->perm, h_jad->perm, n*sizeof(int));
}

/*-------------------------------------------*/
void cuda_malloc_dia(int nd, int strd, 
                     dia_t *d_dia) {
  d_dia->diags = 
  (double*) cuda_malloc(nd*strd*sizeof(double));
  d_dia->ioff = 
  (int*) cuda_malloc(nd*sizeof(int));
}

/*--------------------------------------------*/
void copy_dia_h2d(dia_t *h_dia, dia_t *d_dia) {
  int nd = h_dia->ndiags;
  int strd = h_dia->stride;
  d_dia->n = h_dia->n;
  d_dia->nnz = h_dia->nnz;
  d_dia->ndiags = nd;
  d_dia->stride = strd;
  memcpyh2d(d_dia->diags, h_dia->diags, 
            nd*strd*sizeof(double));
  memcpyh2d(d_dia->ioff, h_dia->ioff, nd*sizeof(int)); 
}


/*-----------------*/
void Free(void *p) {
  if (p) free(p);
}

/*---------------------------*/
void cuda_free(void *p) {
  if (p == NULL) return;
  CUDA_SAFE_CALL(cudaFree(p));
}

/*------------------------------*/
void cuda_free_host(void *p) {
  if (p == NULL) return;
  CUDA_SAFE_CALL(cudaFreeHost(p));
}

/*---------------------------*/
void free_coo(coo_t *coo) {
  if (coo == NULL) return;
  free(coo->ir);
  free(coo->jc);
  free(coo->val);
  free(coo);
}

/*------------------------*/
void free_csr(csr_t *csr) {
  if (csr == NULL) return;
  free(csr->a);
  free(csr->ia);
  free(csr->ja);
  free(csr);
}

/*------------------------------*/
void cuda_free_csr(csr_t *d_csr) {
  if (d_csr == NULL) return;
  cuda_free(d_csr->ia);
  cuda_free(d_csr->ja);
  cuda_free(d_csr->a);
  free(d_csr);
}

/*------------------------*/
void free_jad(jad_t *jad) {
  if (jad == NULL) return;
  free(jad->ia);
  free(jad->ja);
  free(jad->a);
  free(jad->perm);
  free(jad);
}

/*-------------------------------*/
void cuda_free_jad(jad_t *d_jad) {
  if (d_jad == NULL) return;
  cuda_free(d_jad->ia);
  cuda_free(d_jad->ja);
  cuda_free(d_jad->a);
  cuda_free(d_jad->perm);
  cuda_free(d_jad->w);
  free(d_jad);
}

/*------------------------*/
void free_dia(dia_t *dia) {
  if (dia == NULL) return;
  free(dia->diags);
  free(dia->ioff);
  free(dia);
}

/*------------------------------*/
void cuda_free_dia(dia_t *d_dia) {
  if (d_dia == NULL) return;
  cuda_free(d_dia->diags);
  cuda_free(d_dia->ioff);
  free(d_dia);
}


/*-----------------------------*/
void free_matrix(matrix_t *mat) {
  free_csr(mat->h_csr);
  cuda_free_csr(mat->d_csr);
  free_jad(mat->h_jad);
  cuda_free_jad(mat->d_jad);
  free_dia(mat->h_dia);
  cuda_free_dia(mat->d_dia);
  free(mat);
}


