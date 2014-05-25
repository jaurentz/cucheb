/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>


/* chebpoints */
/* complex single precision */
__global__ void cpoints (int n,const cuComplex *a,const cuComplex *b,cuComplex *pts,int incpts){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	float theta;
	cuComplex two = make_cuFloatComplex(2.0f,0.0f);
	cuComplex sine;

	if(ii < n){
		theta = (float)(M_PI_2)*(2*ii-n+1)/(n-1);
		sine = make_cuFloatComplex(sinf(theta),0.0f);
		pts[ii*incpts] = cuCdivf(cuCaddf(cuCaddf(*b,*a),cuCmulf(sine,cuCsubf(*b,*a))),two);
	}
}
cuchebStatus_t cuchebCpoints (int n,const cuComplex *a,const cuComplex *b,cuComplex *pts,int incpts){

	// check n
	if(n <= 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 1.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(n > MAX_FLOAT_DEG+1){
		fprintf(stderr,"\nIn %s line: %d, n must be <= %d.\n",__FILE__,__LINE__,MAX_FLOAT_DEG+1);
		cuchebExit(-1);
	}
	
	// check incpts
	if(incpts <= 0){
		fprintf(stderr,"\nIn %s line: %d, incpts must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// check a and b
	if(&a[0] == &b[0]){
		fprintf(stderr,"\nIn %s line: %d, a must != b.\n",__FILE__,__LINE__);
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
	
	// launch kernel
	cpoints<<<gridSize,blockSize>>>(n,a,b,pts,incpts);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
