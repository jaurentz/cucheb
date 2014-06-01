#include <cucheb.h>


/* chebcoeffs */
/* double precision */
__global__ void dinput (int n, const double *fvals, int incfvals, cufftDoubleReal *input){
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
__global__ void doutput (int n, const cufftDoubleComplex *output, double *coeffs, int inccfs){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii == 0){
		coeffs[ii*inccfs] = cuCreal(output[n-1-ii])/(double)(n-1)/2.0;
	}
	else if(ii == n-1){
		coeffs[ii*inccfs] = cuCreal(output[n-1-ii])/(double)(n-1)/2.0;
	}
	else if(ii > 0 && ii < n-1){
		coeffs[ii*inccfs] = cuCreal(output[n-1-ii])/(double)(n-1);
	}
}
cuchebStatus_t cuchebDcoeffs (int n, const double *fvals, int incfvals, double *coeffs, int inccfs){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(n > MAX_DOUBLE_DEG+1){
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
		cuchebCheckError(cudaMemcpy(coeffs, fvals, sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// n > 1
	else{	
		// allocate workspace
		cufftDoubleReal *input;
		cuchebCheckError(cudaMalloc(&input, 2*(n-1)*sizeof(cufftDoubleReal)),__FILE__,__LINE__);
		cufftDoubleComplex *output;
		cuchebCheckError(cudaMalloc(&output, n*sizeof(cufftDoubleComplex)),__FILE__,__LINE__);
	
		// set blockSize and gridsize
		dim3 blockSize, gridSize;
		cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
		// launch fill input kernel
		dinput<<<gridSize,blockSize>>>(n,fvals,incfvals,input);
	
		// check for kernel error
		cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
		// initialize cufft
		cufftHandle cufftHand;
		cuchebCheckError(cufftPlan1d(&cufftHand, 2*(n-1), CUFFT_D2Z, 1),__FILE__,__LINE__);
	
		// execute plan
		cuchebCheckError(cufftExecD2Z(cufftHand,input,output),__FILE__,__LINE__);
	
		// launch extract output kernel
		doutput<<<gridSize,blockSize>>>(n,output,coeffs,inccfs);
	
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
