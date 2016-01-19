#include <cuchebdependencies.h>

/* header file for cucheblanczos data type */
#ifndef __cucheblanczos_h__ 
#define __cucheblanczos_h__

/* convergence tolerance */
#ifdef DOUBLE_TOL
#undef DOUBLE_TOL
#define DOUBLE_TOL (double)pow(2.0,-52)
#else
#define DOUBLE_TOL (double)pow(2.0,-52)
#endif

/* maximum block size */
#ifdef MAX_BLOCK_SIZE
#undef MAX_BLOCK_SIZE
#define MAX_BLOCK_SIZE 3
#else
#define MAX_BLOCK_SIZE 3
#endif

/* maximum step size */
#ifdef MAX_STEP_SIZE
#undef MAX_STEP_SIZE
#define MAX_STEP_SIZE 6000
#else
#define MAX_STEP_SIZE 6000
#endif

/* maximum number of arnoldi vectors */
#ifdef MAX_NUM_VECS
#undef MAX_NUM_VECS
#define MAX_NUM_VECS 6000
#else
#define MAX_NUM_VECS 6000
#endif

/* maximum number orthogonalization depth */
#ifdef MAX_ORTH_DEPTH
#undef MAX_ORTH_DEPTH
#define MAX_ORTH_DEPTH (MAX_NUM_VECS)
#else
#define MAX_ORTH_DEPTH (MAX_NUM_VECS)
#endif

/* default step size */
#ifdef DEF_STEP_SIZE
#undef DEF_STEP_SIZE
#define DEF_STEP_SIZE 30
#else
#define DEF_STEP_SIZE 30
#endif

/* default number of arnoldi vectors */
#ifdef DEF_NUM_VECS
#undef DEF_NUM_VECS
#define DEF_NUM_VECS 1200
#else
#define DEF_NUM_VECS 1200
#endif

/* cucheblanczos data type */
typedef struct {

  int n;
  int bsize;
  int nblocks;
  int stop;
  int nconv;
  int* index;
  double* evals;
  double* res;
  double* bands;
  double* vecs;
  double* schurvecs;

  double* dtemp;
  double* dvecs;
  double* dschurvecs;
 
} cucheblanczos;

#endif /* __cucheblanczos_h__ */
