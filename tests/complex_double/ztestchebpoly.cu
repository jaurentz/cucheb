#include <cucheb.h>

/* helpers for complex couble precision constructors */
__global__ void zfunkernel(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex *out, int incout);
cuchebStatus_t zfuncaller(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex *out, int incout);

/* drvier */
int main(){

	// compute variables
	int deg;
	cuDoubleComplex a, b;
	double tol;
	
	// ChebPoly 1
	ChebPoly CP(CUCHEB_FIELD_FLOAT_COMPLEX);
	CP.printlong();
	
	// ChebPoly 2
	a = make_cuDoubleComplex(-1.0,0.0);
	b = make_cuDoubleComplex(1.0,0.0);
	deg = pow(2,3);
	ChebPoly CP2(&zfuncaller,&a,&b,deg);
	CP2.printlong();
	
	// ChebPoly 3
	a = make_cuDoubleComplex(-1.0,0.0);
	b = make_cuDoubleComplex(1.0,0.0);
	tol = 1e-7;
	ChebPoly CP3(&zfuncaller,&a,&b,&tol);
	CP3.printlong();
	
	// ChebPoly 4
	a = make_cuDoubleComplex(-1.0,0.0);
	b = make_cuDoubleComplex(1.0,0.0);
	ChebPoly CP4(&zfuncaller,&a,&b);
	CP4.printlong();

	// return	
	return 0;
}

/* helpers for complex double precision constructors */
/* kernel to call device function */
__global__ void zfunkernel(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex *out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		out[ii*incout] = cuCsub(cuCmul(in[ii*incin],in[ii*incin]),in[ii*incin]);
	}
}
/* subroutine to call dfunkernel */
cuchebStatus_t zfuncaller(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex *out, int incout){
	
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
	zfunkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
