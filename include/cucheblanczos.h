#include <cuchebdependencies.h>

/* header file for cucheblanczos data type */
#ifndef __cucheblanczos_h__ 
#define __cucheblanczos_h__

/* maximum number of restarts */
#ifdef MAX_RESTARTS
#undef MAX_RESTARTS
#define MAX_RESTARTS 2
#else
#define MAX_RESTARTS 2
#endif

/* maximum number of arnoldi vectors */
#ifdef MAX_ARNOLDI_VECS
#undef MAX_ARNOLDI_VECS
#define MAX_ARNOLDI_VECS 100
#else
#define MAX_ARNOLDI_VECS 100
#endif

/* maximum number of computed eigenvalues */
#ifdef MAX_NUM_EIGS    
#undef MAX_NUM_EIGS    
#define MAX_NUM_EIGS 100
#else
#define MAX_NUM_EIGS 100
#endif

/* convergence tolerance */
#ifdef DOUBLE_TOL
#undef DOUBLE_TOL
#define DOUBLE_TOL (double)pow(2.0,-52)
#else
#define DOUBLE_TOL (double)pow(2.0,-52)
#endif

/* cucheblanczos data type */
typedef struct {

  int n;
  int nvecs;
  int* index;
  double* diag;
  double* sdiag;
  double* schurvecs;

  double* dtemp;
  double* dvecs;
  double* dschurvecs;
 
} cucheblanczos;

#endif /* __cucheblanczos_h__ */
