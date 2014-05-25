#ifndef __cucheberror_h__ /* __cucheberror_h__ */
#define __cucheberror_h__

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cublas.h>
#include <cublas_v2.h>
#include <curand.h>
#include <cufft.h>
#include <cusparse_v2.h>
#include <math.h>
#include <cuComplex.h>
#include <float.h>

/* c to c++ linker 
#ifdef __cplusplus
extern "C" {
#endif */

/* cucheb exit */
void cuchebExit(int ii);

/* cucheb status type returns */
typedef enum cuchebStatus_t {
    CUCHEB_STATUS_SUCCESS         = 0,
    CUCHEB_STATUS_CUDA_FAILED     = 1,
    CUCHEB_STATUS_CURAND_FAILED   = 2,
    CUCHEB_STATUS_CUFFT_FAILED    = 3,
    CUCHEB_STATUS_CUBLAS_FAILED   = 4,
    CUCHEB_STATUS_CUSPARSE_FAILED = 5
} cuchebStatus_t;

/* error strings */
static const char* cuchebGetErrorString(cuchebStatus_t err);
static const char* cublasGetErrorString(cublasStatus_t err);
static const char* curandGetErrorString(curandStatus_t err);  
static const char* cufftGetErrorString(cufftResult err); 
static const char* cusparseGetErrorString(cusparseStatus_t err); 

/* cucheb error handlers */
void cuchebCheckError(void* err,char* file,int line);            /* host pointers */
void cuchebCheckError(cuchebStatus_t err,char* file,int line);   /* cucheb */
void cuchebCheckError(cudaError_t err,char* file,int line);      /* cuda */
void cuchebCheckError(cublasStatus_t err,char* file,int line);   /* cublas */
void cuchebCheckError(curandStatus_t err,char* file,int line);   /* curand*/
void cuchebCheckError(cufftResult err,char* file,int line);      /* cufft */
void cuchebCheckError(cusparseStatus_t err,char* file,int line); /* cusparse */

/* c to c++ linker 
#ifdef __cplusplus
}
#endif */

#endif /* __cucheberror_h__ */
