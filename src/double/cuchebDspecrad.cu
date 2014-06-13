#include <cucheb.h>

/* computes upper bound on spectral radius interval */
cuchebStatus_t cuchebDspecrad(int n, cuchebOpMult OPMULT, void *USERDATA, double *specrad){

	// check n
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// lanczos parameters
	int runlength = min(60,n-1);

	// allocate memory for Lanzcos
	double *vecs, *diags, *sdiags;
	cuchebCheckError(cudaMalloc(&vecs,(runlength+1)*n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&diags,runlength*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&sdiags,runlength*sizeof(double)),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// initialize curand
	curandGenerator_t curand_gen;
	cuchebCheckError(curandCreateGenerator(&curand_gen, CURAND_RNG_PSEUDO_DEFAULT),__FILE__,__LINE__);
	cuchebCheckError(curandSetPseudoRandomGeneratorSeed(curand_gen,time(NULL)),__FILE__,__LINE__);
	
	// random starting vector
	double temp;
	cuchebCheckError(curandGenerateNormalDouble(curand_gen,vecs,n,0.0,1.0),__FILE__,__LINE__);
//	cuchebCheckError(cuchebDinit(n,vecs,1,1.0),__FILE__,__LINE__);
	cuchebCheckError(cublasDnrm2(cublas_handle,n,vecs,1,&temp),__FILE__,__LINE__);
	temp = 1.0/temp;
	cuchebCheckError(cublasDscal(cublas_handle,n,&temp,vecs,1),__FILE__,__LINE__);
	
	// do Lanzcos run
	cuchebCheckError(cuchebDlanczos(n,OPMULT,USERDATA,0,runlength,vecs,diags,sdiags),__FILE__,__LINE__);

	// initialize gpu eigenvectors
	double *eigvecs;
	cuchebCheckError(cudaMalloc(&eigvecs,runlength*runlength*sizeof(double)),__FILE__,__LINE__);

	// initialize cpu diags
	double *h_diags;
	cuchebCheckError((void*)(h_diags = (double*)malloc(runlength*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(h_diags,diags,runlength*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

	// initialize cpu sdiags
	double *h_sdiags;
	cuchebCheckError((void*)(h_sdiags = (double*)malloc((runlength-1)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(h_sdiags,sdiags,(runlength-1)*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

	// initialize cpu eigenvectors
	double *h_eigvecs;
	cuchebCheckError((void*)(h_eigvecs = (double*)malloc(runlength*runlength*sizeof(double))),__FILE__,__LINE__);
	for(int ii=0;ii<runlength;ii++){
		for(int jj=0;jj<runlength;jj++){
			if(ii == jj){
				h_eigvecs[ii+jj*runlength] = 1.0;
			}
			else{
				h_eigvecs[ii+jj*runlength] = 0.0;
			}
		}
	}

	// call lapacke
	lapack_int lpint = 0;
	lpint = LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I',runlength,h_diags,h_sdiags,h_eigvecs,runlength);
	if(lpint |= 0){
		fprintf(stderr,"\nLapacke error: %d occured in %s at line: %d\n\n",lpint,__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// copy results back to gpu
	cuchebCheckError(cudaMemcpy(diags,h_diags,runlength*sizeof(double),cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(eigvecs,h_eigvecs,runlength*runlength*sizeof(double),cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute ritz vectors
	double *ritz_min, *ritz_max;
	cuchebCheckError(cudaMalloc(&ritz_min,n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ritz_max,n*sizeof(double)),__FILE__,__LINE__);
	temp = 1.0;
	double temp2 = 0.0;
	cuchebCheckError(cublasDgemv(cublas_handle,CUBLAS_OP_N,n,runlength,&temp,vecs,n,&eigvecs[0],1,&temp2,ritz_min,1),__FILE__,__LINE__);
	cuchebCheckError(cublasDgemv(cublas_handle,CUBLAS_OP_N,n,runlength,&temp,vecs,n,&eigvecs[runlength*(runlength-1)],1,&temp2,ritz_max,1),__FILE__,__LINE__);
	cuchebCheckError(cublasDnrm2(cublas_handle,n,ritz_min,1,&temp),__FILE__,__LINE__);
	temp = 1.0/temp;
	cuchebCheckError(cublasDscal(cublas_handle,n,&temp,ritz_min,1),__FILE__,__LINE__);
	cuchebCheckError(cublasDnrm2(cublas_handle,n,ritz_max,1,&temp),__FILE__,__LINE__);
	temp = 1.0/temp;
	cuchebCheckError(cublasDscal(cublas_handle,n,&temp,ritz_max,1),__FILE__,__LINE__);
	
	// compute ritz values
	double r1, r2;
	OPMULT(ritz_min,&vecs[0],USERDATA);
	cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[0],1,ritz_min,1,&r1),__FILE__,__LINE__);
	OPMULT(ritz_max,&vecs[n],USERDATA);
	cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[n],1,ritz_max,1,&r2),__FILE__,__LINE__);
	
	// compute residuals
	temp = -r1;
	cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,ritz_min,1,&vecs[0],1),__FILE__,__LINE__);
	temp = -r2;
	cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,ritz_max,1,&vecs[n],1),__FILE__,__LINE__);
	double res1, res2;
	cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[0],1,&res1),__FILE__,__LINE__);
	cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[n],1,&res2),__FILE__,__LINE__);
	
	// set specrad if converged
	double tol = 1e-2;
	int max_restarts = 10;
	*specrad = max(abs(r1),abs(r2));
	if(abs(*specrad) == 0){
		fprintf(stderr,"\nIn %s line: %d, specrad must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(abs(r1) < abs(r2) && res2 < tol*(*specrad)){
		*specrad = abs(r2*1.1);
	}
	else if(abs(r2) < abs(r1) && res1 < tol*(*specrad)){
		*specrad = abs(r1*1.1);
	}
	// restart if not converged
	else{
		if(abs(r1) > abs(r2)){
			cuchebCheckError(cudaMemcpy(ritz_max,ritz_min,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		}
		for(int ii=0;ii<max_restarts;ii++){
			// starting vector
			cuchebCheckError(cuchebDinit(n*(runlength+1),vecs,1,0.0),__FILE__,__LINE__);
			cuchebCheckError(cudaMemcpy(&vecs[0],ritz_max,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	
			// do Lanzcos run
			cuchebCheckError(cuchebDlanczos(n,OPMULT,USERDATA,0,runlength,vecs,diags,sdiags),__FILE__,__LINE__);

			// initialize cpu diags
			cuchebCheckError(cudaMemcpy(h_diags,diags,runlength*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

			// initialize cpu sdiags
			cuchebCheckError(cudaMemcpy(h_sdiags,sdiags,(runlength-1)*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

			// initialize cpu eigenvectors
			cuchebCheckError(cudaMemcpy(h_eigvecs,eigvecs,runlength*runlength*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			for(int kk=0;kk<runlength;kk++){
				for(int jj=0;jj<runlength;jj++){
					if(kk == jj){
						h_eigvecs[kk+jj*runlength] = 1.0;
					}
					else{
						h_eigvecs[kk+jj*runlength] = 0.0;
					}
				}
			}

			// call lapacke
			lapack_int lpint;
			lpint = LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I',runlength,h_diags,h_sdiags,h_eigvecs,runlength);
			if(lpint |= 0){
				fprintf(stderr,"\nLapacke error: %d occured in %s at line: %d\n\n",lpint,__FILE__,__LINE__);
				cuchebExit(-1);
			}

			// copy results back to gpu
			cuchebCheckError(cudaMemcpy(diags,h_diags,runlength*sizeof(double),cudaMemcpyHostToDevice),__FILE__,__LINE__);
			cuchebCheckError(cudaMemcpy(eigvecs,h_eigvecs,runlength*runlength*sizeof(double),cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
			// compute ritz vectors
			temp = 1.0;
			temp2 = 0.0;
			cuchebCheckError(cublasDgemv(cublas_handle,CUBLAS_OP_N,n,runlength,&temp,vecs,n,&eigvecs[runlength*(runlength-1)],1,&temp2,ritz_max,1),__FILE__,__LINE__);
			cuchebCheckError(cublasDnrm2(cublas_handle,n,ritz_max,1,&temp),__FILE__,__LINE__);
			temp = 1.0/temp;
			cuchebCheckError(cublasDscal(cublas_handle,n,&temp,ritz_max,1),__FILE__,__LINE__);
	
			// compute ritz values
			OPMULT(ritz_max,&vecs[n],USERDATA);
			cuchebCheckError(cublasDdot(cublas_handle,n,&vecs[n],1,ritz_max,1,&r2),__FILE__,__LINE__);
	
			// compute residuals
			temp = -r2;
			cuchebCheckError(cublasDaxpy(cublas_handle,n,&temp,ritz_max,1,&vecs[n],1),__FILE__,__LINE__);
			cuchebCheckError(cublasDnrm2(cublas_handle,n,&vecs[n],1,&res2),__FILE__,__LINE__);

			// check convergence
			if(res2 < tol*abs(r2)){
				*specrad = abs(r2*1.1);
			}
			
			// check iterations
			if(ii == max_restarts-1){
				fprintf(stderr,"\nIn %s line: %d, restarts maxed out.\n",__FILE__,__LINE__);
				cuchebExit(-1);
			}
		}	
	}
	
			
	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);
	
	// shutdown curand
	cuchebCheckError(curandDestroyGenerator(curand_gen),__FILE__,__LINE__);

	// free memory
	cuchebCheckError(cudaFree(vecs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(diags),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(sdiags),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(eigvecs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ritz_min),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ritz_max),__FILE__,__LINE__);

	free(h_eigvecs);
	free(h_diags);
	free(h_sdiags);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}
