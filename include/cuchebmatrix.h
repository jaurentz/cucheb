/* header file for cuchebmatrix data type */
#ifndef __cuchebmatrix_h__ 
#define __cuchebmatrix_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
using namespace std;

#include <mmio.h>
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuchebpoly.h>

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

/* instantiate cuchebmatrix object */
int cuchebmatrix_init(const string& mtxfile, cuchebmatrix* ccm);

/* destroy cuchebmatrix object */
int cuchebmatrix_destroy(cuchebmatrix* ccm);

/* print cuchebmatrix object */
int cuchebmatrix_print(cuchebmatrix* ccm);

/* longprint cuchebmatrix object */
int cuchebmatrix_printlong(cuchebmatrix* ccm);

/* gpuprint cuchebmatrix object */
int cuchebmatrix_gpuprint(cuchebmatrix* ccm);

/* routine for sorting entries */
int cuchebmatrix_sort(cuchebmatrix* ccm);

/* routine for converting to csr format */
int cuchebmatrix_csr(cuchebmatrix* ccm);

/* routine for mv multiply on GPU */
int cuchebmatrix_mv(cuchebmatrix* ccm, double* alpha, double* x, double* beta,
                    double* y);

/* routine for poly mv multiply on GPU */
int cuchebmatrix_polymv(cuchebmatrix* ccm, cuchebpoly* ccp, double* alpha, double* x, 
                        double* beta, double* y);

/* routine for estimating spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm);

#endif /* __cuchebmatrix_h__ */
