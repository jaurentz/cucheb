#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <datatype.h>
#include <protos.h>
#include <cuda.h>

void compare_mv(matrix_t *mat);

/*-------------------------------*/
int main() {
/*-------------------------------*/
/*   driver program of GPU-sol   */
/*-------------------------------*/
  char finput[] = "input";
  options_t *opts;
  matrix_t  *mat;
  coo_t     *h_coo;
  int       i, n, nnz;
/*----------------------- init */
  if (gpusol_init()) {
    printf("GPU-Solv Init Error\n");
    return 1;
  }
/*------------ read input file */
  Calloc(opts, 1, options_t);
  read_input(finput, opts);
/*----------- first step: ----------*/
/*-- alloc 'mat' & 'prec' structure */
/*----------------------------------*/
  Calloc(mat, 1, matrix_t);
/* load matrix,stored in COO struct */
  Calloc(h_coo, 1, coo_t);
  read_coo_MM(h_coo, opts);
  n = h_coo->n;  nnz = h_coo->nnz;
/*-------- convert to CSR in 'mat' */
  mat->n = n;  mat->nnz = nnz;
  Calloc(mat->h_csr, 1, csr_t);
  COO2CSR(h_coo, mat->h_csr);
  setup_matrix(mat, opts);

  compare_mv(mat);

/*------ done, free mem */
  free_opts(opts);
  free_matrix(mat);                                                                                                                                                                                      
  free_coo(h_coo);
/*------finalize & check error */

}






