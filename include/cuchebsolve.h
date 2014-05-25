#ifndef __cuchebsolve_h__ /* __cuchebsolve_h__ */
#define __cuchebsolve_h__

#include <cuchebop.h>

cuchebStatus_t cuchebDsolve(int n, cuchebOpMult OPMULT, void* USERDATA, double *x, double *b);

#endif /* __cuchebsolve_h__ */
