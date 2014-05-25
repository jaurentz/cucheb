/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>


/* chebcoeffs */
/* complex double precision */
__global__ void zinput (int n, const cuDoubleComplex *fvals, int incfvals, cufftDoubleComplex *input){
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
__global__ void zoutput (int n, const cufftDoubleComplex *output, cuDoubleComplex *coeffs, int inccfs){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	cuDoubleComplex scl;

	if(ii == 0){
		scl = make_cuDoubleComplex((double)(n-1)*2.0,0.0);
		coeffs[ii*inccfs] = cuCdiv(output[n-1-ii],scl);
	}
	else if(ii == n-1){
		scl = make_cuDoubleComplex((double)(n-1)*2.0,0.0);
		coeffs[ii*inccfs] = cuCdiv(output[n-1-ii],scl);
	}
	else if(ii > 0 && ii < n-1){
		scl = make_cuDoubleComplex((double)(n-1),0.0);
		coeffs[ii*inccfs] = cuCdiv(output[n-1-ii],scl);
	}
}
cuchebStatus_t cuchebZcoeffs (int n, const cuDoubleComplex *fvals, int incfvals, cuDoubleComplex *coeffs, int inccfs){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	else if(n > MAX_DOUBLE_DEG+1){
		fprintf(stderr,"\nIn %s line: %d, n must be <= %d.\n",__FILE__,__LINE__,MAX_DOUBLE_DEG+1);
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
		cuchebCheckError(cudaMemcpy(coeffs, fvals, sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// n > 1
	else{
		// allocate workspace
		cufftDoubleComplex *input;
		cuchebCheckError(cudaMalloc(&input, 2*(n-1)*sizeof(cufftDoubleComplex)),__FILE__,__LINE__);
	
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
		zinput<<<gridSize,blockSize>>>(n,fvals,incfvals,input);
	
		// check for kernel error
		cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
		// initialize cufft
		cufftHandle cufftHand;
		cuchebCheckError(cufftPlan1d(&cufftHand, 2*(n-1), CUFFT_Z2Z, 1),__FILE__,__LINE__);
	
		// execute plan
		cuchebCheckError(cufftExecZ2Z(cufftHand,input,input,CUFFT_FORWARD),__FILE__,__LINE__);
	
		// launch extract output kernel
		zoutput<<<gridSize,blockSize>>>(n,input,coeffs,inccfs);
	
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
