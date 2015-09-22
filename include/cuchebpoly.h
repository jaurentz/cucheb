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

/* maximum construction degree */
#ifdef DOUBLE_DEG
#undef DOUBLE_DEG
#define DOUBLE_DEG 8192
#else
#define DOUBLE_DEG 8192
#endif

/* construction tolerance */
#ifdef DOUBLE_EPS
#undef DOUBLE_EPS
#define DOUBLE_EPS (double)pow(2.0,-8)
#else
#define DOUBLE_EPS (double)pow(2.0,-8)
#endif

/* maximum filter degree */
#ifdef MAX_DOUBLE_DEG
#undef MAX_DOUBLE_DEG
#define MAX_DOUBLE_DEG 1024
#else
#define MAX_DOUBLE_DEG 1024
#endif

/* cuchebpoly data type */
typedef struct {

  int degree;
  double a;
  double b;
  double points[2*DOUBLE_DEG];
  double coeffs[DOUBLE_DEG+1];

  cufftHandle cuffthandle;
  cufftDoubleReal *dinput;
  cufftDoubleComplex *doutput;

} cuchebpoly;

#endif /* __cuchebpoly_h__ */
