#include <cucheb.h>

// prototypes for sample functions to interpolate
/* complex double precision */
__host__ __device__ void ztestfun(const cuDoubleComplex *x, cuDoubleComplex *y);
__global__ void ztestkernel(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex* out, int incout);
cuchebStatus_t ztestcaller(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex* out, int incout);

// driver
int main(){

	// compute variables
	/* complex double */
	int zn = pow(2,16)+1;
	int sizeZ = sizeof(cuDoubleComplex);
	cuDoubleComplex *zpts, *za, *zb, *zcfs, *zfvs, *dzpts, *dza, *dzb, *dzcfs, *dzfvs;
	
	// allocate host memory
	/* complex double */
	cuchebCheckError((void*)(zpts = (cuDoubleComplex*)malloc(zn*sizeZ)),__FILE__,__LINE__);
	cuchebCheckError((void*)(zcfs = (cuDoubleComplex*)malloc(zn*sizeZ)),__FILE__,__LINE__);
	cuchebCheckError((void*)(zfvs = (cuDoubleComplex*)malloc(zn*sizeZ)),__FILE__,__LINE__);
	cuchebCheckError((void*)(za = (cuDoubleComplex*)malloc(sizeZ)),__FILE__,__LINE__);
	cuchebCheckError((void*)(zb = (cuDoubleComplex*)malloc(sizeZ)),__FILE__,__LINE__);
	
	// allocate device memory
	/* complex double */
	cuchebCheckError(cudaMalloc(&dzpts, zn*sizeZ),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dzcfs, zn*sizeZ),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dzfvs, zn*sizeZ),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dza, sizeZ),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dzb, sizeZ),__FILE__,__LINE__);
	
	// set host pointers
	/* complex double */
	*za = make_cuDoubleComplex(-1.0,-1.0);
	*zb = make_cuDoubleComplex(1.0,3.0);
	
	// copy host memory to device memory
	/* complex double */
	cuchebCheckError(cudaMemcpy(dza, za, sizeZ, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dzb, zb, sizeZ, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* complex double */
	cuchebCheckError(cuchebZpoints(zn,dza,dzb,dzpts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex double */
	cuchebCheckError(cudaMemcpy(zpts, dzpts, zn*sizeZ, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute funvals
	/* complex double */
	cuchebCheckError(ztestcaller(zn, dzpts, 1, dzfvs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex double */
	cuchebCheckError(cudaMemcpy(zfvs, dzfvs, zn*sizeZ, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute chebcoeffs
	/* complex double */
	cuchebCheckError(cuchebZcoeffs(zn, dzfvs, 1, dzcfs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex double */
	cuchebCheckError(cudaMemcpy(zcfs, dzcfs, zn*sizeZ, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* complex double */
	printf("complex double precision\n");
	for(int ii=0;ii<8;ii++){printf("zfvs[%d] = (%+1.15e,%+1.15e), zcfs[%d] = (%+1.15e,%+1.15e)\n",
		ii,cuCreal(zfvs[ii]),cuCimag(zfvs[ii]),ii,cuCreal(zcfs[ii]),cuCimag(zcfs[ii]));}
	printf("\n");

	// free device memory
	/* complex double */
	cuchebCheckError(cudaFree(dzpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dzcfs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dzfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dza),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dzb),__FILE__,__LINE__);

	// free host memory
	/* complex double */
	free(zpts);
	free(zcfs);
	free(zfvs);
	free(za);
	free(zb);
	
	return 0;
}

// sample functions to interpolate
/* complex double precision */
__host__ __device__ void ztestfun(const cuDoubleComplex *x, cuDoubleComplex *y){
	//cuDoubleComplex temp;
	//cuDoubleComplex I = make_cuDoubleComplex(1.0,1.0);
	//temp = cuCmul(*x,I);
	*y = *x;
}
__global__ void ztestkernel(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex* out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		ztestfun(&in[ii*incin],&out[ii*incout]);
	}
}
cuchebStatus_t ztestcaller(int n, const cuDoubleComplex *in, int incin, cuDoubleComplex* out, int incout){
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
	ztestkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
