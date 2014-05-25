#include <cucheb.h>

// prototypes for sample functions to interpolate
/* complex single precision */
__host__ __device__ void ctestfun(const cuComplex *x, cuComplex *y);
__global__ void ctestkernel(int n, const cuComplex *in, int incin, cuComplex* out, int incout);
cuchebStatus_t ctestcaller(int n, const cuComplex *in, int incin, cuComplex* out, int incout);

// driver
int main(){

	// compute variables
	/* complex single */
	int cn = pow(2,3)+1;
	int sizeC = sizeof(cuComplex);
	cuComplex *cpts, *ca, *cb, *ccfs, *cfvs, *dcpts, *dca, *dcb, *dccfs, *dcfvs;
	
	// allocate host memory
	/* complex single */
	cuchebCheckError((void*)(cpts = (cuComplex*)malloc(cn*sizeC)),__FILE__,__LINE__);
	cuchebCheckError((void*)(ccfs = (cuComplex*)malloc(cn*sizeC)),__FILE__,__LINE__);
	cuchebCheckError((void*)(cfvs = (cuComplex*)malloc(cn*sizeC)),__FILE__,__LINE__);
	cuchebCheckError((void*)(ca = (cuComplex*)malloc(sizeC)),__FILE__,__LINE__);
	cuchebCheckError((void*)(cb = (cuComplex*)malloc(sizeC)),__FILE__,__LINE__);
	
	// allocate device memory
	/* complex single */
	cuchebCheckError(cudaMalloc(&dcpts, cn*sizeC),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dccfs, cn*sizeC),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dcfvs, cn*sizeC),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dca, sizeC),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dcb, sizeC),__FILE__,__LINE__);
	
	// set host pointers
	/* complex single */
	*ca = make_cuFloatComplex(-1.0f,0.0f);
	*cb = make_cuFloatComplex(1.0f,0.0f);
	
	// copy host memory to device memory
	/* complex single */
	cuchebCheckError(cudaMemcpy(dca, ca, sizeC, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dcb, cb, sizeC, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* complex single */
	cuchebCheckError(cuchebCpoints(cn,dca,dcb,dcpts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex single */
	cuchebCheckError(cudaMemcpy(cpts, dcpts, cn*sizeC, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute funvals
	/* complex single */
	cuchebCheckError(ctestcaller(cn, dcpts, 1, dcfvs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex single */
	cuchebCheckError(cudaMemcpy(cfvs, dcfvs, cn*sizeC, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// compute chebcoeffs
	/* complex single */
	cuchebCheckError(cuchebCcoeffs(cn, dcfvs, 1, dccfs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex single */
	cuchebCheckError(cudaMemcpy(ccfs, dccfs, cn*sizeC, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* complex single */
	printf("complex single precision\n");
	for(int ii=0;ii<cn;ii++){printf("cpts[%d] = (%+e,%+e), cfvs[%d] = (%+e,%+e), ccfs[%d] = (%+e,%+e)\n",
		ii,cuCrealf(cpts[ii]),cuCimagf(cpts[ii]),ii,cuCrealf(cfvs[ii]),cuCimagf(cfvs[ii]),ii,cuCrealf(ccfs[ii]),cuCimagf(ccfs[ii]));}
	printf("\n");

	// free device memory
	/* complex single */
	cuchebCheckError(cudaFree(dcpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dccfs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dcfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dca),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dcb),__FILE__,__LINE__);

	// free host memory
	/* complex single */
	free(cpts);
	free(ccfs);
	free(cfvs);
	free(ca);
	free(cb);
	
	return 0;
}

// sample functions to interpolate
/* complex single precision */
__host__ __device__ void ctestfun(const cuComplex *x, cuComplex *y){
	cuComplex temp;
	cuComplex I = make_cuFloatComplex(1.0f,1.0f);
	temp = cuCmulf(*x,I);
	*y = cuCmulf(cuCmulf(temp,temp),temp);
}
__global__ void ctestkernel(int n, const cuComplex *in, int incin, cuComplex* out, int incout){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < n){
		ctestfun(&in[ii*incin],&out[ii*incout]);
	}
}
cuchebStatus_t ctestcaller(int n, const cuComplex *in, int incin, cuComplex* out, int incout){
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
	ctestkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
