#include <cuchebdependencies.h>

/* header file for cuchebmatrix data type */
#ifndef __cuchebmatrix_h__ 
#define __cuchebmatrix_h__

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
