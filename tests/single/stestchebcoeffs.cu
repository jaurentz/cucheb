#include <cucheb.h>

// prototypes for sample functions to interpolate
/* single precision */
__host__ __device__ void stestfun(const float *x, float *y);
__global__ void stestkernel(int n, const float *in, int incin, float* out, int incout);
cuchebStatus_t stestcaller(int n, const float *in, int incin, float* out, int incout);

// driver
int main(){

	// compute variables
	/* single */
	int sn = pow(2,3)+1;
	int sizeS = sizeof(float);
	float *spts, *sa, *sb, *scfs, *sfvs, *dspts, *dsa, *dsb, *dscfs, *dsfvs;
	
	// allocate host memory
	/* single */
	cuchebCheckError((void*)(spts = (float*)malloc(sn*sizeS)),__FILE__,__LINE__);
	cuchebCheckError((void*)(scfs = (float*)malloc(sn*sizeS)),__FILE__,__LINE__);
	cuchebCheckError((void*)(sfvs = (float*)malloc(sn*sizeS)),__FILE__,__LINE__);
	cuchebCheckError((void*)(sa = (float*)malloc(sizeS)),__FILE__,__LINE__);
	cuchebCheckError((void*)(sb = (float*)malloc(sizeS)),__FILE__,__LINE__);
	
	// allocate device memory
	/* single */
	cuchebCheckError(cudaMalloc(&dspts, sn*sizeS),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dscfs, sn*sizeS),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dsfvs, sn*sizeS),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dsa, sizeS),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dsb, sizeS),__FILE__,__LINE__);
	
	// set host pointers
	/* single */
	*sa = -1.0f;
	*sb = 1.0f;
	
	// copy host memory to device memory
	/* single */
	cuchebCheckError(cudaMemcpy(dsa, sa, sizeS, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dsb, sb, sizeS, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* single */
	cuchebCheckError(cuchebSpoints(sn,dsa,dsb,dspts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* single */
	cuchebCheckError(cudaMemcpy(spts, dspts, sn*sizeS, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute funvals
	/* single */
	cuchebCheckError(stestcaller(sn, dspts, 1, dsfvs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* single */
	cuchebCheckError(cudaMemcpy(sfvs, dsfvs, sn*sizeS, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute chebcoeffs
	/* single */
	cuchebCheckError(cuchebScoeffs(sn, dsfvs, 1, dscfs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* single */
	cuchebCheckError(cudaMemcpy(scfs, dscfs, sn*sizeS, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* single */
	printf("single precision\n");
	for(int ii=0;ii<sn;ii++){printf("spts[%d] = %+e, sfvs[%d] = %+e, scfs[%d] = %+e\n",ii,spts[ii],ii,sfvs[ii],ii,scfs[ii]);}
	printf("\n");

	// free device memory
	/* single */
	cuchebCheckError(cudaFree(dspts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dscfs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dsfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dsa),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dsb),__FILE__,__LINE__);

	// free host memory
	/* single */
	free(spts);
	free(scfs);
	free(sfvs);
	free(sa);
	free(sb);
	
	return 0;
}

// sample functions to interpolate
/* single precision */
__host__ __device__ void stestfun(const float *x, float *y){
	*y = expf(*x);
}
__global__ void stestkernel(int n, const float *in, int incin, float* out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		stestfun(&in[ii*incin],&out[ii*incout]);
	}
}
cuchebStatus_t stestcaller(int n, const float *in, int incin, float* out, int incout){
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
	stestkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
