#include <cucheb.h>

// prototypes for sample functions to interpolate
/* double precision */
__host__ __device__ void dtestfun(const double *x, double *y);
__global__ void dtestkernel(int n, const double *in, int incin, double *out, int incout);
cuchebStatus_t dtestcaller(int n, const double *in, int incin, double *out, int incout);

// driver
int main(){

	// compute variables
	/* double */
	int dn = pow(2,3)+1;
	int sizeD = sizeof(double);
	double *dpts, *da, *db, *dcfs, *dfvs, *ddpts, *dda, *ddb, *ddcfs, *ddfvs;
	
	// allocate host memory
	/* double */
	cuchebCheckError((void*)(dpts = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(dcfs = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(dfvs = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(da = (double*)malloc(sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(db = (double*)malloc(sizeD)),__FILE__,__LINE__);
	
	// allocate device memory
	/* double */
	cuchebCheckError(cudaMalloc(&ddpts, dn*sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ddcfs, dn*sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ddfvs, dn*sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dda, sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ddb, sizeD),__FILE__,__LINE__);
	
	// set host pointers
	/* double */
	*da = -1.0;
	*db = 1.0;
	
	// copy host memory to device memory
	/* double */
	cuchebCheckError(cudaMemcpy(dda, da, sizeD, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(ddb, db, sizeD, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* double */
	cuchebCheckError(cuchebDpoints(dn,dda,ddb,ddpts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* double */
	cuchebCheckError(cudaMemcpy(dpts, ddpts, dn*sizeD, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute funvals
	/* double */
	cuchebCheckError(dtestcaller(dn, ddpts, 1, ddfvs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* double */
	cuchebCheckError(cudaMemcpy(dfvs, ddfvs, dn*sizeD, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute chebcoeffs
	/* double */
	cuchebCheckError(cuchebDcoeffs(dn, ddfvs, 1, ddcfs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* double */
	cuchebCheckError(cudaMemcpy(dcfs, ddcfs, dn*sizeD, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* double */
	printf("double precision\n");
	for(int ii=0;ii<dn;ii++){printf("dpts[%d] = %+1.15e, dfvs[%d] = %+1.15e, dcfs[%d] = %+1.15e\n",ii,dpts[ii],ii,dfvs[ii],ii,dcfs[ii]);}
	printf("\n");

	// free device memory
	/* double */
	cuchebCheckError(cudaFree(ddpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ddcfs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ddfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dda),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ddb),__FILE__,__LINE__);

	// free host memory
	/* double */
	free(dpts);
	free(dcfs);
	free(dfvs);
	free(da);
	free(db);
	
	return 0;
}

// sample functions to interpolate
/* double precision */
__host__ __device__ void dtestfun(const double *x, double *y){
	*y = sin(*x);
}
__global__ void dtestkernel(int n, const double *in, int incin, double* out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		dtestfun(&in[ii*incin],&out[ii*incout]);
	}
}
cuchebStatus_t dtestcaller(int n, const double *in, int incin, double* out, int incout){
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
	dtestkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
