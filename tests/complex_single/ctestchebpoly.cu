#include <cucheb.h>

/* helpers for cuComplex precision constructors */
__global__ void cfunkernel(int n, const cuComplex *in, int incin, cuComplex *out, int incout);
cuchebStatus_t cfuncaller(int n, const cuComplex *in, int incin, cuComplex *out, int incout);

/* drvier */
int main(){

	// compute variables
	int deg;
	cuComplex a, b;
	float tol;
	
	// ChebPoly 1
	ChebPoly CP(CUCHEB_FIELD_FLOAT_COMPLEX);
	CP.printlong();
	
	// ChebPoly 2
	a = make_cuFloatComplex(-1.0f,0.0f);
	b = make_cuFloatComplex(1.0f,0.0f);
	deg = pow(2,3);
	ChebPoly CP2(&cfuncaller,&a,&b,deg);
	CP2.printlong();
	
	// ChebPoly 2
	a = make_cuFloatComplex(-1.0f,0.0f);
	b = make_cuFloatComplex(1.0f,0.0f);
	tol = 1e-5;
	ChebPoly CP3(&cfuncaller,&a,&b,&tol);
	CP3.print();
	
	// ChebPoly 4
	a = make_cuFloatComplex(-1.0f,0.0f);
	b = make_cuFloatComplex(1.0f,0.0f);
	ChebPoly CP4(&cfuncaller,&a,&b);
	CP4.print();

	// return	
	return 0;
}

/* helpers for cuComplex precision constructors */
/* kernel to call device function */
__global__ void cfunkernel(int n, const cuComplex *in, int incin, cuComplex *out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		out[ii*incout] = make_cuFloatComplex(123.0f*cuCabsf(in[ii*incin]),0.0f);
	}
}
/* subroutine to call dfunkernel */
cuchebStatus_t cfuncaller(int n, const cuComplex *in, int incin, cuComplex *out, int incout){
	
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
	cfunkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
