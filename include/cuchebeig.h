#ifndef __cuchebeig_h__ /* __cuchebeig_h__ */
#define __cuchebeig_h__

#include <cuchebop.h>

cuchebStatus_t cuchebDlanczos(int n,cuchebOpMult OPMULT,void *USERDATA,int start,int runlength,double *vecs,double *diags,double *sdiags);
cuchebStatus_t cuchebDlanczos(ChebOp *CP, int start, int runlength, double *vecs, double *diags, double *sdiags);
cuchebStatus_t cuchebDspecrad(int n, cuchebOpMult OPMULT, void *USERDATA, double *specrad);
cuchebStatus_t cuchebDrestart(int n,int runlength,int neigs,int *nconv,double *vecs,double *diags,double *sdiags,double *ritzvecs,double tol);

#endif /* __cuchebeig_h__ */
