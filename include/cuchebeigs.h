#ifndef __cuchebeigs_h__
#define __cuchebeigs_h__

#include <cuchebop.h>

struct cuchebLanczosHandle{
	int n;
	int numeigs;
	int runlength;
	int restarts;
	double tol;
	int numconv;
	int numrestarts;
	int nummatvecs;
};

cuchebStatus_t cuchebDeigs(cuchebLanczosHandle* LH, cuchebOpMult OPMULT, void* USERDATA, double *eigvecs);
cuchebStatus_t cuchebDeigs(cuchebLanczosHandle* LH, ChebOp * CO, double *eigvecs);

cuchebStatus_t cuchebDlanczos(int n,cuchebOpMult OPMULT,void *USERDATA,int start,int runlength,double *vecs,double *diags,double *sdiags);
cuchebStatus_t cuchebDlanczos(ChebOp *CP, int start, int runlength, double *vecs, double *diags, double *sdiags);
cuchebStatus_t cuchebDspecrad(int n, cuchebOpMult OPMULT, void *USERDATA, double *specrad);
cuchebStatus_t cuchebDrestart(int n,int runlength,int neigs,int *nconv,double *vecs,double *diags,double *sdiags,double *ritzvecs,double tol);

#endif
