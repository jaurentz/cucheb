#include <cucheb.h>

/* approximation to 1/z */
__global__ void invDkernel(int n, const double *in, int incin, double *out, int incout, double *userdata);
cuchebStatus_t invDcaller(int n, const double *in, int incin, double *out, int incout, void *userdata);

/* linear system solver */
cuchebStatus_t cuchebDsolve(int n, cuchebOpMult OPMULT, void *USERDATA, double *x, double *b){

	// check n
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// compute spect interval
	double A,B;
//	cuchebCheckError(cuchebDsolve(n,OPMULT,USERDATA,&A,&B),__FILE__,__LINE__);
	
	// compute approximate inverse 
	double tol = 1e-10;
	double temp;
	double *tau;
	cuchebCheckError(cudaMalloc(&tau,sizeof(double)),__FILE__,__LINE__);
	temp = -log(tol)/A;
	cuchebCheckError(cuchebDinit(1,tau,1,temp),__FILE__,__LINE__);
	ChebPoly CP(&invDcaller,&A,&B,tau,&tol);
	//CP.print();
	
	// set chebop
	ChebOp CO(n,OPMULT,USERDATA,&CP);
	CO.print();
	
	// solve Ax = b
	cuchebCheckError(CO.Mult(b,x),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// allocate memory for Lanzcos
	double *vecs;
	cuchebCheckError(cudaMalloc(&vecs,3*n*sizeof(double)),__FILE__,__LINE__);
	
	// compute residual
	cuchebCheckError(cudaMemcpy(&vecs[0],b,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	double nrmb;
	cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[0],1,&nrmb),__FILE__,__LINE__);
	OPMULT(x,&vecs[n],USERDATA);
	temp = -1.0;
	cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[n],1,&vecs[0],1),__FILE__,__LINE__);
	double nrmres;
	cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[0],1,&nrmres),__FILE__,__LINE__);
	printf("res = %+e\n",nrmres/nrmb);
	
	// iterative refinement
	int max_refinements = 1;
	for(int ii=0;ii<max_refinements;ii++){
	
		// correction
		cuchebCheckError(CO.Mult(&vecs[0],&vecs[2*n]),__FILE__,__LINE__);
		temp = 1.0;
		cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[2*n],1,x,1),__FILE__,__LINE__);
		
		// new residual
		cuchebCheckError(cudaMemcpy(&vecs[0],b,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		OPMULT(x,&vecs[n],USERDATA);
		temp = -1.0;
		cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[n],1,&vecs[0],1),__FILE__,__LINE__);
		cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[0],1,&nrmres),__FILE__,__LINE__);
		printf("res = %+e\n",nrmres/nrmb);
		printf("\n");
	}

	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);

	// free memory
	cuchebCheckError(cudaFree(tau),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(vecs),__FILE__,__LINE__);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}

/* approximation to 1/z */
__global__ void invDkernel(int n, const double *in, int incin, double *out, int incout, double *userdata){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	if(ii < n){out[ii*incout] = (1.0 - exp(-(*userdata)*in[ii*incin]))/in[ii*incin];}
}
cuchebStatus_t invDcaller(int n, const double *in, int incin, double *out, int incout, void *userdata){
	
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
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// launch fill input kernel
	invDkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout, (double*)userdata);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
