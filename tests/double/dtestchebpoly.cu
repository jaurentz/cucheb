#include <cucheb.h>

/* helpers for double precision constructors */
__global__ void dfunkernel(int n, const double *in, int incin, double *out, int incout);
cuchebStatus_t dfuncaller(int n, const double *in, int incin, double *out, int incout);

/* drvier */
int main(){

	// compute variables
	int deg;
	double a, b;
	double tol;
	
	// ChebPoly 1
	ChebPoly CP(CUCHEB_FIELD_DOUBLE);
	CP.printlong();
	
	// ChebPoly 2
	a = -1.0;
	b = 1.0;
	deg = pow(2,4);
	ChebPoly CP2(&dfuncaller,&a,&b,deg);
	CP2.printlong();
	
	// ChebPoly 3
	a = -1.0;
	b = 1.0;
	tol = 1e-3;
	ChebPoly CP3(&dfuncaller,&a,&b,&tol);
	CP3.print();
	
	// ChebPoly 4
	a = -1.0;
	b = 1.0;
	ChebPoly CP4(&dfuncaller,&a,&b);
	CP4.print();

	// return	
	return 0;
}

/* helpers for double precision constructors */
/* kernel to call device function */
__global__ void dfunkernel(int n, const double *in, int incin, double *out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		//out[ii*incout] = exp(sin(10.0*M_PI_2*in[ii*incin]));
		out[ii*incout] = abs(in[ii*incin]);
	}
}
/* subroutine to call dfunkernel */
cuchebStatus_t dfuncaller(int n, const double *in, int incin, double *out, int incout){
	
	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check incin
	if(incin <= 0){
		fprintf(stderr,"\nIn %s line: %d, incin must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check incout
	if(incout <= 0){
		fprintf(stderr,"\nIn %s line: %d, incout must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// query device
	int dev;
	cudaDeviceProp prop;
	cuchebCheckError(cudaGetDevice(&dev),__FILE__,__LINE__);
	cuchebCheckError(cudaGetDeviceProperties(&prop,dev),__FILE__,__LINE__);
	
	// set blockSize
	int blockSize;
	blockSize = prop.maxThreadsPerBlock;
	
	// set gridSize
	int gridSize;
	gridSize = (int)ceil((double)n/blockSize);
	
	// launch fill input kernel
	dfunkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
