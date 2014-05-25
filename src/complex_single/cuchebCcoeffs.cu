/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>


/* chebcoeffs */
/* complex single precision */
__global__ void cinput (int n, const cuComplex *fvals, int incfvals, cufftComplex *input){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
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
__global__ void coutput (int n, const cufftComplex *output, cuComplex *coeffs, int inccfs){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	cuComplex scl;

	if(ii == 0){
		scl = make_cuFloatComplex((float)(n-1)*2.0f,0.0f);
		coeffs[ii*inccfs] = cuCdivf(output[n-1-ii],scl);
	}
	else if(ii == n-1){
		scl = make_cuFloatComplex((float)(n-1)*2.0f,0.0f);
		coeffs[ii*inccfs] = cuCdivf(output[n-1-ii],scl);
	}
	else if(ii > 0 && ii < n-1){
		scl = make_cuFloatComplex((float)(n-1),0.0f);
		coeffs[ii*inccfs] = cuCdivf(output[n-1-ii],scl);
	}
}
cuchebStatus_t cuchebCcoeffs (int n, const cuComplex *fvals, int incfvals, cuComplex *coeffs, int inccfs){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	else if(n > MAX_FLOAT_DEG+1){
		fprintf(stderr,"\nIn %s line: %d, n must be <= %d.\n",__FILE__,__LINE__,MAX_FLOAT_DEG+1);
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
		cuchebCheckError(cudaMemcpy(coeffs, fvals, sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// n > 1
	else{
		// allocate workspace
		cufftComplex *input;
		cuchebCheckError(cudaMalloc(&input, 2*(n-1)*sizeof(cufftComplex)),__FILE__,__LINE__);
	
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
		cinput<<<gridSize,blockSize>>>(n,fvals,incfvals,input);
	
		// check for kernel error
		cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
		// initialize cufft
		cufftHandle cufftHand;
		cuchebCheckError(cufftPlan1d(&cufftHand, 2*(n-1), CUFFT_C2C, 1),__FILE__,__LINE__);
	
		// execute plan
		cuchebCheckError(cufftExecC2C(cufftHand,input,input,CUFFT_FORWARD),__FILE__,__LINE__);
	
		// launch extract output kernel
		coutput<<<gridSize,blockSize>>>(n,input,coeffs,inccfs);
	
		// check for kernel error
		cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
		// free cufft
		cuchebCheckError(cufftDestroy(cufftHand),__FILE__,__LINE__);
	
		// free workspace
		cuchebCheckError(cudaFree(input),__FILE__,__LINE__);
	}
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
