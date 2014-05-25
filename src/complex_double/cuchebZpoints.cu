/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>


/* chebpoints */
/* complex double precision */
__global__ void zpoints (int n,const cuDoubleComplex *a,const cuDoubleComplex *b,cuDoubleComplex *pts,int incpts){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	double theta;
	cuDoubleComplex two = make_cuDoubleComplex(2.0,0.0);
	cuDoubleComplex sine;

	if(ii < n){
		theta = (double)(M_PI_2)*(2*ii-n+1)/(n-1);
		sine = make_cuDoubleComplex(sin(theta),0.0);
		pts[ii*incpts] = cuCdiv(cuCadd(cuCadd(*b,*a),cuCmul(sine,cuCsub(*b,*a))),two);
	}
}
cuchebStatus_t cuchebZpoints (int n,const cuDoubleComplex *a,const cuDoubleComplex *b,cuDoubleComplex *pts,int incpts){

	// check n
	if(n <= 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 1.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(n > MAX_DOUBLE_DEG+1){
		fprintf(stderr,"\nIn %s line: %d, n must be <= %d.\n",__FILE__,__LINE__,MAX_DOUBLE_DEG+1);
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
	zpoints<<<gridSize,blockSize>>>(n,a,b,pts,incpts);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
