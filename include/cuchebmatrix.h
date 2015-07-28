/* header file for cuchebmatrix data type */
#ifndef __cuchebmatrix_h__ 
#define __cuchebmatrix_h__

#include <stdio.h>
#include <stdlib.h>
#include <mmio.h>
#include <cuda.h>
#include <cusparse.h>

/* cuchebmatrix data type */
typedef struct {

  MM_typecode matcode;
  int m;
  int n;
  int nnz;
  int* rowinds;
  int* colinds;
  double* vals;

} cuchebmatrix;

/* instantiate cuchebmatrix object */
int cuchebmatrix_init(char* mtxfile, cuchebmatrix* ccm);

/* destroy cuchebmatrix object */
int cuchebmatrix_destroy(cuchebmatrix* ccm);

/* print cuchebmatrix object */
int cuchebmatrix_print(cuchebmatrix* ccm);

/* longprint cuchebmatrix object */
int cuchebmatrix_printlong(cuchebmatrix* ccm);

/* routine for sorting entries using GPU */
int cuchebmatrix_sort(cuchebmatrix* ccm);

#endif /* __cuchebmatrix_h__ */
