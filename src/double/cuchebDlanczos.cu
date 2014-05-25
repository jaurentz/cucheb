#include <cucheb.h>

cuchebStatus_t cuchebDlanczos(int n, cuchebOpMult OPMULT, void *USERDATA, int start, int runlength, double *vecs, double *diags, double *sdiags){ 
	
	// check start
	if(start < 0){
		fprintf(stderr,"\nIn %s line: %d, start must be => 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check runlength
	if(runlength < 1){
		fprintf(stderr,"\nIn %s line: %d, runlength must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check n
	if(n <= (start+runlength)){
		fprintf(stderr,"\nIn %s line: %d, n must be > start+runlength.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
		
	// lanczos run
	double temp;
	for(int ii=start;ii<runlength;ii++){	
		// compute A*v
		OPMULT((void*)&vecs[ii*n],(void*)&vecs[(ii+1)*n],USERDATA);

		// orthogonalize
		for(int jj=0;jj<(ii+1);jj++){
			// dot prod
			cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[(ii+1)*n],1,&vecs[(ii-jj)*n],1,&temp),__FILE__,__LINE__);

			// store diag
			if(jj == 0){cuchebCheckError(cublasSetVector(1,sizeof(double),&temp,1,&diags[ii],1),__FILE__,__LINE__);}

			// subtract
			temp = -temp;
			cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[(ii-jj)*n],1,&vecs[(ii+1)*n],1),__FILE__,__LINE__);
		}

		// reorthogonalize
		for(int jj=0;jj<(ii+1);jj++){

			// dot prod
			cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[(ii+1)*n],1,&vecs[(ii-jj)*n],1,&temp),__FILE__,__LINE__);

			// subract
			temp = -temp;
			cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[(ii-jj)*n],1,&vecs[(ii+1)*n],1),__FILE__,__LINE__);
		}

		// compute norm
		cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[(ii+1)*n],1,&temp),__FILE__,__LINE__);

		// store sdiag
		cuchebCheckError(cublasSetVector(1,sizeof(double),&temp,1,&sdiags[ii],1),__FILE__,__LINE__);

		// normalize
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,&vecs[(ii+1)*n],1),__FILE__,__LINE__);
	}
	
	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);	

	// return
	return CUCHEB_STATUS_SUCCESS;
}

cuchebStatus_t cuchebDlanczos(ChebOp *CP, int start, int runlength, double *vecs, double *diags, double *sdiags){ 
	
	// check start
	if(start < 0){
		fprintf(stderr,"\nIn %s line: %d, start must be => 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check runlength
	if(runlength < 1){
		fprintf(stderr,"\nIn %s line: %d, runlength must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check n
	int n = CP->getN();
	if(n <= (start+runlength)){
		fprintf(stderr,"\nIn %s line: %d, n must be > start+runlength.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
		
	// lanczos run
	double temp;
	for(int ii=start;ii<runlength;ii++){	
		// compute A*v
		cuchebCheckError(CP->Mult(&vecs[ii*n],&vecs[(ii+1)*n]),__FILE__,__LINE__);

		// orthogonalize
		for(int jj=0;jj<(ii+1);jj++){
			// dot prod
			cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[(ii+1)*n],1,&vecs[(ii-jj)*n],1,&temp),__FILE__,__LINE__);

			// store diag
			if(jj == 0){cuchebCheckError(cublasSetVector(1,sizeof(double),&temp,1,&diags[ii],1),__FILE__,__LINE__);}

			// subtract
			temp = -temp;
			cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[(ii-jj)*n],1,&vecs[(ii+1)*n],1),__FILE__,__LINE__);
		}

		// reorthogonalize
		for(int jj=0;jj<(ii+1);jj++){

			// dot prod
			cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[(ii+1)*n],1,&vecs[(ii-jj)*n],1,&temp),__FILE__,__LINE__);

			// subract
			temp = -temp;
			cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,&vecs[(ii-jj)*n],1,&vecs[(ii+1)*n],1),__FILE__,__LINE__);
		}

		// compute norm
		cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[(ii+1)*n],1,&temp),__FILE__,__LINE__);

		// store sdiag
		cuchebCheckError(cublasSetVector(1,sizeof(double),&temp,1,&sdiags[ii],1),__FILE__,__LINE__);

		// normalize
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,&vecs[(ii+1)*n],1),__FILE__,__LINE__);
	}
	
	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);	

	// return
	return CUCHEB_STATUS_SUCCESS;
}
