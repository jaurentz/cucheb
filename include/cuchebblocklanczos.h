#include <cuchebdependencies.h>

/* header file for cuchebblocklanczos data type */
#ifndef __cuchebblocklanczos_h__ 
#define __cuchebblocklanczos_h__

/* maximum number of restarts */
#ifdef MAX_BLOCK_SIZE
#undef MAX_BLOCK_SIZE
#define MAX_BLOCK_SIZE 3
#else
#define MAX_BLOCK_SIZE 3
#endif

/* maximum number of arnoldi vectors */
#ifdef MAX_NUM_BLOCKS
#undef MAX_NUM_BLOCKS
#define MAX_NUM_BLOCKS 100
#else
#define MAX_NUM_BLOCKS 100
#endif

/* cuchebblocklanczos data type */
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
 
} cuchebblocklanczos;

#endif /* __cuchebblocklanczos_h__ */
