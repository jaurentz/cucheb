#include <cucheb.h>

__global__ void dinit(int n,double *x,int incx,double val){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

	if(ii < n){
		x[ii*incx] = val;
	}
}

cuchebStatus_t cuchebDinit(int n,double *x,int incx,double val){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check incx
	if(incx <= 0){
		fprintf(stderr,"\nIn %s line: %d, incx must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);

	// call kernel
	dinit<<<gridSize,blockSize>>>(n,x,incx,val);
	//dinit<<<1,n>>>(n,x,incx,val);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return
	return CUCHEB_STATUS_SUCCESS;
}
