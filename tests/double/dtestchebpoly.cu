#include <cucheb.h>

/* helpers for double precision constructors */
cuchebStatus_t dfuncaller(int n, const double *in, int incin, double *out, int incout, void* userdata){

	double *tau = (double*)userdata;
	
	for(int ii=0;ii<n;ii++){
		out[ii*incout] = sin((*tau)*in[ii*incin]);
	}
	
	return CUCHEB_STATUS_SUCCESS;

}

/* driver */
int main(){

	// compute variables
	int deg;
	double a, b;
	double tol;
	double scl = 4.0*M_PI_2;
	void* userdata = &scl;
	
	// ChebPoly 1
	ChebPoly CP(CUCHEB_FIELD_DOUBLE);
	CP.printlong();
	
	// ChebPoly 2
	a = -1.0;
	b = 1.0;
	deg = pow(2,4);
	ChebPoly CP2(&dfuncaller,userdata,&a,&b,deg);
	CP2.printlong();
	
	// ChebPoly 3
	a = -1.0;
	b = 1.0;
	tol = 1e-5;
	ChebPoly CP3(&dfuncaller,userdata,&a,&b,&tol);
	CP3.printlong();
	
	// ChebPoly 4
	a = -1.0;
	b = 1.0;
	ChebPoly CP4(&dfuncaller,userdata,&a,&b);
	CP4.printlong();

	// return	
	return 0;
}
