#include <cuchebdependencies.h>

/* header file for cuchebpoly data type */
#ifndef __cuchebpoly_h__ 
#define __cuchebpoly_h__

/* double precision pi */
#ifdef DOUBLE_PI
#undef DOUBLE_PI
#define DOUBLE_PI 3.141592653589793238
#else
#define DOUBLE_PI 3.141592653589793238
#endif

/* construction tolerance */
#ifdef DOUBLE_EPS
#undef DOUBLE_EPS
#define DOUBLE_EPS (double)pow(2.0,-6)
#else
#define DOUBLE_EPS (double)pow(2.0,-6)
#endif

/* maximum filter degree */
#ifdef MAX_DOUBLE_DEG
#undef MAX_DOUBLE_DEG
#define MAX_DOUBLE_DEG 8192
#else
#define MAX_DOUBLE_DEG 8192
#endif

/* cuchebpoly data type */
typedef struct {

  int degree;
  double a;
  double b;
  double points[2*MAX_DOUBLE_DEG];
  double coeffs[MAX_DOUBLE_DEG+1];

  cufftHandle cuffthandle;
  cufftDoubleReal *dinput;
  cufftDoubleComplex *doutput;

} cuchebpoly;

#endif /* __cuchebpoly_h__ */
