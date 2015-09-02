#include <cuchebdependencies.h>

/* header file for cucheblanczos data type */
#ifndef __cucheblanczos_h__ 
#define __cucheblanczos_h__

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

#endif /* __cucheblanczos_h__ */
