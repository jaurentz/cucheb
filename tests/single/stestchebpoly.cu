#include <cucheb.h>

/* helpers for single precision constructors */
__global__ void sfunkernel(int n, const float *in, int incin, float* out, int incout);
cuchebStatus_t sfuncaller(int n, const float *in, int incin, float* out, int incout, void* userdata);

/* driver */
int main(){

	// compute variables
	int deg;
	float a, b;
	float tol;
	
	// ChebPoly 1
	ChebPoly CP(CUCHEB_FIELD_FLOAT);
	CP.printlong();
	
	// ChebPoly 2
	a = -1.0f;
	b = 1.0f;
	deg = pow(2,4);
	ChebPoly CP2(&sfuncaller,&a,&b,deg);
	CP2.printlong();
	
	// ChebPoly 3
	a = -1.0f;
	b = 1.0f;
	tol = 1e-3;
	ChebPoly CP3(&sfuncaller,&a,&b,&tol);
	CP3.printlong();
	
	// ChebPoly 4
	a = -1.0f;
	b = 1.0f;
	ChebPoly CP4(&sfuncaller,&a,&b);
	CP4.printlong();

	// return	
	return 0;
}

/* helpers for single precision constructors */
/* kernel to call devince function */
__global__ void sfunkernel(int n, const float *in, int incin, float* out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		out[ii*incout] = expf(sinf(4.0f*M_PI_2*in[ii*incin]));
		//out[ii*incout] = 101.0f*in[ii*incin];
	}
}
/* subroutine to call sfunkernel */
cuchebStatus_t sfuncaller(int n, const float *in, int incin, float* out, int incout, void* userdata){
	
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
	sfunkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
