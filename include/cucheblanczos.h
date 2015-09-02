/* header file for cucheblanczos data type */
#ifndef __cucheblanczos_h__ 
#define __cucheblanczos_h__

#include <lapacke.h>
#include <cuchebmatrix.h>

/* cucheblanczos data type */
typedef struct {

  int n;
  int nvecs;
  double* diag;
  double* sdiag;
  double* schurvecs;

  double* dtemp;
  double* dvecs;
  double* dschurvecs;
 
} cucheblanczos;

/* instantiate cucheblanczos object */
int cucheblanczos_init(int nvecs, cuchebmatrix* ccm, cucheblanczos* ccl);

/* destroy cucheblanczos object */
int cucheblanczos_destroy(cucheblanczos* ccl);

/* print cucheblanczos object */
int cucheblanczos_print(cucheblanczos* ccl);

/* set cucheblanczos starting vector */
int cucheblanczos_startvec(cucheblanczos* ccl);

/* arnoldi run using cuchebmatrix */
int cucheblanczos_arnoldi(cuchebmatrix* ccm, cucheblanczos* ccl);

/* compute ritz values and vectors */
int cucheblanczos_eig(cuchebmatrix* ccm, cucheblanczos* ccl);

#endif /* __cucheblanczos_h__ */
