#include <cucheb.h>

/* chebcoeffs */
/* single precision */
__global__ void sinput (int n, const float *fvals, int incfvals, cufftReal *input){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	int N = 2*(n-1);

	if(ii == 0){
		input[ii] = fvals[(n-1)*incfvals];
	}
	else if(ii == n-1){
		input[ii] = fvals[0];
	}
	else if(ii > 0 && ii < n-1){
		input[ii] = fvals[(n-1-ii)*incfvals];
		input[N-ii] = fvals[(n-1-ii)*incfvals];
	}
}
__global__ void soutput (int n, const cufftComplex *output, float *coeffs, int inccfs){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

	if(ii == 0){
		coeffs[ii*inccfs] = cuCrealf(output[n-1-ii])/(float)(n-1)/2.0f;
	}
	else if(ii == n-1){
		coeffs[ii*inccfs] = cuCrealf(output[n-1-ii])/(float)(n-1)/2.0f;
	}
	else if(ii > 0 && ii < n-1){
		coeffs[ii*inccfs] = cuCrealf(output[n-1-ii])/(float)(n-1);
	}
}
cuchebStatus_t cuchebScoeffs (int n, const float *fvals, int incfvals, float *coeffs, int inccfs){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(n > (MAX_FLOAT_DEG+1)){
		fprintf(stderr,"\nIn %s line: %d, n must be <= %d.\n",__FILE__,__LINE__,(MAX_FLOAT_DEG+1));
		cuchebExit(-1);
	}
	
	// check incfvals
	if(incfvals <= 0){
		fprintf(stderr,"\nIn %s line: %d, incfvals must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check inccfs
	if(inccfs <= 0){
		fprintf(stderr,"\nIn %s line: %d, inccfs must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// n = 1
	if(n == 1){
		cuchebCheckError(cudaMemcpy(coeffs, fvals, sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// n > 1
	else{
		// allocate workspace
		cufftReal *input;
		cuchebCheckError(cudaMalloc(&input, 2*(n-1)*sizeof(cufftReal)),__FILE__,__LINE__);
		cufftComplex *output;
		cuchebCheckError(cudaMalloc(&output, n*sizeof(cufftComplex)),__FILE__,__LINE__);
	
		// set blockSize and gridsize
		dim3 blockSize, gridSize;
		cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
		// launch fill input kernel
		sinput<<<gridSize,blockSize>>>(n,fvals,incfvals,input);
	
		// check for kernel error
		cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
		// initialize cufft
		cufftHandle cufftHand;
		cuchebCheckError(cufftPlan1d(&cufftHand, 2*(n-1), CUFFT_R2C, 1),__FILE__,__LINE__);
	
		// execute plan
		cuchebCheckError(cufftExecR2C(cufftHand,input,output),__FILE__,__LINE__);
	
		// launch extract output kernel
		soutput<<<gridSize,blockSize>>>(n,output,coeffs,inccfs);
	
		// check for kernel error
		cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
		// free cufft
		cuchebCheckError(cufftDestroy(cufftHand),__FILE__,__LINE__);
	
		// free workspace
		cuchebCheckError(cudaFree(input),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(output),__FILE__,__LINE__);
	}
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
