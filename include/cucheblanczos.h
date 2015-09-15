#include <cuchebdependencies.h>

/* header file for cucheblanczos data type */
#ifndef __cucheblanczos_h__ 
#define __cucheblanczos_h__

/* convergence tolerance */
#ifdef DOUBLE_TOL
#undef DOUBLE_TOL
#define DOUBLE_TOL (double)pow(2.0,-40)
#else
#define DOUBLE_TOL (double)pow(2.0,-40)
#endif

/* maximum number of restarts */
#ifdef MAX_RESTARTS
#undef MAX_RESTARTS
#define MAX_RESTARTS 2
#else
#define MAX_RESTARTS 2
#endif

/* maximum number of computed eigenvalues */
#ifdef MAX_NUM_EIGS    
#undef MAX_NUM_EIGS    
#define MAX_NUM_EIGS 100
#else
#define MAX_NUM_EIGS 100
#endif

/* maximum number of restarts */
#ifdef MAX_BLOCK_SIZE
#undef MAX_BLOCK_SIZE
#define MAX_BLOCK_SIZE 10
#else
#define MAX_BLOCK_SIZE 10
#endif

/* maximum number of arnoldi vectors */
#ifdef MAX_NUM_BLOCKS
#undef MAX_NUM_BLOCKS
#define MAX_NUM_BLOCKS 100
#else
#define MAX_NUM_BLOCKS 100
#endif

/* maximum number of arnoldi vectors */
#ifdef MAX_ORTH_DEPTH
#undef MAX_ORTH_DEPTH
#define MAX_ORTH_DEPTH 100
#else
#define MAX_ORTH_DEPTH 100
#endif

/* cucheblanczos data type */
typedef struct {

  int n;
  int bsize;
  int nblocks;
  int* index;
  double* evals;
  double* res;
  double* bands;
  double* schurvecs;

  double* dtemp;
  double* dvecs;
  double* dschurvecs;
 
} cucheblanczos;

#endif /* __cucheblanczos_h__ */
