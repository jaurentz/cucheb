#include <cucheb.h>

/* chebpoints */
/* double precision */
__global__ void dpoints (int n,const double *a,const double *b,double *pts,int incpts){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	double theta;

	if(ii < n){
		theta = (double)(M_PI_2)*(2*ii-n+1)/(n-1);
		pts[ii*incpts] = (*b + *a)/2.0 + sin(theta)*(*b - *a)/2.0;
	}
}
cuchebStatus_t cuchebDpoints (int n,const double *a,const double *b,double *pts,int incpts){

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
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// launch kernel
	dpoints<<<gridSize,blockSize>>>(n,a,b,pts,incpts);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
