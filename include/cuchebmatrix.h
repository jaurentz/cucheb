#include <cuchebdependencies.h>

/* header file for cuchebmatrix data type */
#ifndef __cuchebmatrix_h__ 
#define __cuchebmatrix_h__
/*

  cuchebmatrix

  This file defines the cuchebmatrix object. This object is needed for storing
  and using sparse matrices. When an instance of a cuchebmatrix is initialized
  the entries are first put into the proper order on the CPU. Once this is done
  a second copy of the matrix is created on the GPU. The second copy on the GPU
  is necessary in order to avoid making time consuming memory swaps between the
  CPU and GPU.

  CPU variables:

    matcode - matcode used in Matrix Market files
    m       - number of rows 
    n       - number of columns
    nnz     - number of nonzero entries
    a,b     - upper and lower bounds on the eigenvalues
    rowinds - pointer to array of row indices
    colinds - pointer to array of column indices
    vals    - pointer to array of nonzero entries

  GPU variables:

    cublashandle   - handle needed to use CUBLAS subroutines
    cusparsehandle - handle needed to use CUSPARSE subroutines
    matdescr       - used to describe matrix storage on GPU
    drowinds       - pointer to array of row indices stored on GPU
    dcolinds       - pointer to array of column indices stored on GPU
    dvals          - pointer to array of nonzero entries stored on GPU
    dtemp          - pointer to array used for swap space on GPU

*/

/* cuchebmatrix data type */
typedef struct {

  MM_typecode matcode;
  int m;
  int n;
  int nnz;
  double a, b;
  int* rowinds;
  int* colinds;
  double* vals;

  cublasHandle_t cublashandle;
  cusparseHandle_t cusparsehandle;
  cusparseMatDescr_t matdescr;
  int* drowinds;
  int* dcolinds;
  double* dvals;
  double* dtemp;
 
} cuchebmatrix;

#endif /* __cuchebmatrix_h__ */
