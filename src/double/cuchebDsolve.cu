#include <cucheb.h>

/* approximation to 1/z */
cuchebStatus_t invDcaller(int n, const double *in, int incin, double *out, int incout, void* userdata){

	double *tau = (double*)userdata;
	
	for(int ii=0;ii<n;ii++){
		if(in[ii*incin] != 0){
			out[ii*incout] = (1.0-exp(-(*tau)*in[ii*incin]))/in[ii*incin];
		}
		else{
			out[ii*incout] = (*tau);
		}
	}
	
	return CUCHEB_STATUS_SUCCESS;

}

/* linear system solver */
cuchebStatus_t cuchebDsolve(int n, cuchebOpMult OPMULT, void *USERDATA, double *x, double *b){

	// check n
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// compute spect interval
	double A,B;
	A = 0.0;
	B = 1.0;
	
	// compute approximate inverse 
	double tol = 1e-15;
	double tau = 1e6;
	ChebPoly CP(&invDcaller,&tau,&A,&B,&tol);
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
	
	// allocate memory for residual
	double *vecs;
	cuchebCheckError(cudaMalloc(&vecs,3*n*sizeof(double)),__FILE__,__LINE__);
	
	// compute norm b
	cuchebCheckError(cudaMemcpy(&vecs[0],b,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	double nrmb;
	cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[0],1,&nrmb),__FILE__,__LINE__);
	
	// compute Ax
	OPMULT(x,&vecs[n],USERDATA);
	
	// compute residual
	double temp = -1.0;
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
	cuchebCheckError(cudaFree(vecs),__FILE__,__LINE__);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}
