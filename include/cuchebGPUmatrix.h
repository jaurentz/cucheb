/* header file for cuchebGPUmatrix data type */
#ifndef __cuchebGPUmatrix_h__ 
#define __cuchebGPUmatrix_h__

#include <stdio.h>
#include <stdlib.h>
#include <mmio.h>
#include <cuda.h>
#include <cusparse.h>
#include <cuchebmatrix.h>

/* cuchebGPUmatrix data type */
typedef struct {

  cusparseMatDescr_t* matdescr;
  int m;
  int n;
  int nnz;
  int* d_rowinds;
  int* d_colinds;
  double* d_vals;

} cuchebGPUmatrix;

/* instantiate cuchebGPUmatrix object */
int cuchebGPUmatrix_init(cuchebmatrix* ccm, cuchebGPUmatrix* ccGm);

/* destroy cuchebGPUmatrix object */
int cuchebGPUmatrix_destroy(cuchebGPUmatrix* ccGm);

#endif /* __cuchebGPUmatrix_h__ */
