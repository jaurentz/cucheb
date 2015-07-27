#include "gpusollib.h"

/**/
void COO2CSR(coo_t *coo, csr_t *csr) {

   int n, nnz;
   // Allocate CSR
   n = coo->n;
   nnz = coo->nnz;
   malloc_csr(n, nnz, csr);
   // COO to CSR
   coo2csr(n, nnz, coo->val, coo->ir, coo->jc, csr->a, csr->ja, csr->ia);
}

void coo2csr(int n, int nnz, double *coo_vals, int *coo_rows, int *coo_cols, double *acsr, int *ia, int *ja) {
  // The COO vectors hold A in a column-wise manner
  // This implementation is inneficient

  int  i, j, counter1 = -1;
  
  ia[0] = 0;

  for ( i = 0; i < n; i++) {
     for ( j = 0; j < nnz; j++ ) {
       if (coo_rows[j] == i) {
	  counter1++;
          acsr[counter1] = coo_vals[j];
          ja[counter1]   = coo_cols[j];
       }
    }
     ia[i+1] = counter1 + 1;
  }

}


/*------------------------------------------*/
void PadJAD32(jad_t *jad) {
  int i,*oldia,*oldja,njad,nnz2,jadlen;
  double *olda;
/*-----------------------------------*/
  oldia = jad->ia;
  oldja = jad->ja;
  olda  = jad->a;
  njad = jad->njad;
  Malloc(jad->ia, njad+1, int);
  jad->ia[0] = 1;
/*------ pad each jad multiple of 32 */
  for (i=0; i<njad; i++) {
     jadlen = (oldia[i+1]-oldia[i]+31)/32*32;
     jad->ia[i+1] = jad->ia[i] + jadlen;
  }
/*------------------------------*/
  nnz2 = jad->ia[njad]-1;
  jad->nnz = nnz2;
  Calloc(jad->a, nnz2, double);
  Malloc(jad->ja, nnz2, int);
  for (i=0; i<nnz2; i++)
    jad->ja[i] = 1;
/*-------- copy jad */
  for (i=0; i<njad; i++) {
    int p = jad->ia[i] - 1;
    int q = oldia[i] - 1;
    memcpy(&jad->a[p], &olda[q], 
    (oldia[i+1]-oldia[i])*sizeof(double)); 
    memcpy(&jad->ja[p], &oldja[q],
    (oldia[i+1]-oldia[i])*sizeof(int));    
  }
/*----------*/
  free(olda);
  free(oldia);
  free(oldja);
}

/*----------------------------------------*/
void CSR2JAD(csr_t *csr, jad_t *jad) {
/*---------- Allocate JAD */
  jad->n = csr->n;
  jad->nnz = csr->nnz;
  Malloc(jad->ia, csr->n+1, int);
  Malloc(jad->ja, csr->nnz, int);
  Malloc(jad->a, csr->nnz, double);
  Malloc(jad->perm, csr->n, int);
/*---------- CSR -> JAD */
  csrjad_(&csr->n, csr->a, csr->ja, csr->ia, 
          &jad->njad, jad->perm, jad->a, 
          jad->ja, jad->ia);
  assert(jad->njad < BLOCKDIM);
/*--------- pad each jad to multiple of 32 */
  PadJAD32(jad);
}

/*--------------------------------------------*/
int CSR2DIA(csr_t *csr, dia_t *dia) {
  int n = csr->n;
  int ndiags = dia->ndiags = MAXDIAG;
/*------------------ Allocate DIA */
  dia->n = n;
/*------- pad each diag to be multiple of 32 */
  dia->stride = (n+31)/32*32;
  Malloc(dia->diags, dia->stride*ndiags, double);
  Malloc(dia->ioff, ndiags, int);
/*------------ job code = 10: 
   select diags internally, no fill for remainder
   check skit.f for details */
  int job = 10;
/*--------------- work array */
  int *ind = (int*) malloc((2*n-1)*sizeof(int));
  csrdia_(&n, &(dia->ndiags), &job, csr->a, csr->ja,
          csr->ia, &(dia->stride), dia->diags,
          dia->ioff, NULL, NULL, NULL, ind);
  free(ind);
/*-------------------------------*/
  if (dia->ndiags < ndiags) {
    ndiags = dia->ndiags;
    Realloc(dia->diags, dia->stride*ndiags, double);
    Realloc(dia->ioff,ndiags, int);
    dia->nnz = csr->nnz;
    return 1;
  }
/*---- fail to convert */
  return 0;
}

/*-------------------------------------------------*/
void csrcsc(int n, int n2, int job, int ipos, 
            double *a, int *ja, int *ia, 
	    double *ao, int *jao, int *iao) {
  int i,j,k,next;
/*------- compute lengths of rows of A' */
  for (i=1; i<=n2+1; i++)
    iao[i-1] = 0;
    
  for (i=1; i<=n; i++)
    for (k=ia[i-1]; k<=ia[i]-1; k++) {
      j = ja[k-1]+1;
      iao[j-1] ++;
    }
/*---- compute pointers from lengths */  
  iao[0] = ipos;
  for (i=1; i<=n2; i++)
    iao[i] += iao[i-1];
/*---- now do the actual copying */
  for (i=1; i<=n; i++)
    for (k=ia[i-1]; k<=ia[i]-1; k++) {
      j = ja[k-1];
      next = iao[j-1];
      if (job == 1)
        ao[next-1] = a[k-1];
      jao[next-1] = i;
      iao[j-1] = next + 1;
    }
/*---- reshift iao and leave */
  for (i=n2; i>=1; i--)
    iao[i] = iao[i-1];
  iao[0] = ipos;
}

/*-------------------------------------------*
 * Sort each row by increasing column order
 * By double transposition
 *-------------------------------------------*/
void sortrow(csr_t *A) {
  int n, nnz, *ia, *ja;
  double *a;
/*-------------------------------------------*/
  n = A->n;
  ia = A->ia;
  ja = A->ja;
  a = A->a;
  nnz = ia[n] - 1;
  // work array
  double *b;
  Malloc(b, nnz, double);
  int *jb, *ib;
  Malloc(jb, nnz, int);
  Malloc(ib, n+1, int);
  // double transposition
  csrcsc(n, n, 1, 1, a, ja, ia, b, jb, ib);
  csrcsc(n, n, 1, 1, b, jb, ib, a, ja, ia);
  free(b);
  free(jb);
  free(ib);
}

/*-----------------------------------------------*/
void setup_matrix(matrix_t *mat, int choice) {

  int n, nnz, njad, nd, strd, err;
  int csr_flag, jad_flag, dia_flag;
/*-----------------------------------------------*/
  n = mat->n;
  csr_flag = jad_flag = dia_flag = 0;
/*-----------------------------------------------*/
  
  switch ( choice ) {
  case 0: // CSR
    csr_flag = 1;
  break;
  case 1: // JAD
    jad_flag = 1;
    mat->spmv = spmv_jad;
  break;
  case 2: // DIA
    dia_flag = 1;
    mat->spmv = spmv_dia;
  }

/*---------------------------------*/
/* sort rows by increasing col idx */
/*---------------------------------*/
//  sortrow(mat->h_csr);

/*-------------------------------------*/
  if (csr_flag) {
/*----------- copy csr to device */
    Calloc(mat->d_csr, 1, csr_t);
    nnz = mat->h_csr->nnz;
    cuda_malloc_csr(n, nnz, mat->d_csr);
    copy_csr_h2d(mat->h_csr, mat->d_csr);
  }
/*-------------------------------------*/
  if (jad_flag) {
/*------- convert to jad format */
    Calloc(mat->h_jad, 1, jad_t);
    CSR2JAD(mat->h_csr, mat->h_jad);
/*-------- copy jad to device */
    Calloc(mat->d_jad, 1, jad_t);
    nnz = mat->h_jad->nnz; /*NOTE: nnz changed*/
    njad = mat->h_jad->njad;
    cuda_malloc_jad(n, njad, nnz, mat->d_jad);
    copy_jad_h2d(mat->h_jad, mat->d_jad);
  }
/*-------------------------------------*/
  if (dia_flag) {
/*------- convert to dia format */
    Calloc(mat->h_dia, 1, dia_t);
    err = CSR2DIA(mat->h_csr, mat->h_dia);
    assert(err == 1);
/*-------- copy dia to device */
    Calloc(mat->d_dia, 1, dia_t);
    nd = mat->h_dia->ndiags;
    strd = mat->h_dia->stride;
    cuda_malloc_dia(nd, strd, mat->d_dia);
    copy_dia_h2d(mat->h_dia, mat->d_dia);
  }

}

