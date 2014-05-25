#include <cucheb.h>

/* complex cuDoubleComplex precision constructors */
/* fixed degree */
ChebPoly::ChebPoly(cuchebCuDoubleComplexFun fun, cuDoubleComplex *A, cuDoubleComplex *B, void *USERDATA, int Deg){

	// set field
	field = CUCHEB_FIELD_DOUBLE_COMPLEX;

	// check degree
	if(Deg < 0){
		fprintf(stderr,"\nIn %s line: %d, degree must be >= 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(Deg > MAX_DOUBLE_DEG){
		fprintf(stderr,"\nIn %s line: %d, degree must be <= %d.\n",__FILE__,__LINE__,MAX_DOUBLE_DEG);
		cuchebExit(-1);
	}
	
	// set degree
	degree = Deg;
	
	// check a and b
	if(cuCreal(*A) == cuCreal(*B) && cuCimag(*A) == cuCimag(*B)){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// degree 0
	if(degree == 0){
		// compute funvals
		cuDoubleComplex *sfvs;
		cuchebCheckError(cudaMalloc(&sfvs, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError((*fun)(1, (cuDoubleComplex*)a, 1, sfvs, 1, USERDATA),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, sfvs, sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);	
	}
	
	// degree > 0
	else{
		// compute chebpoints
		cuDoubleComplex *spts;
		cuchebCheckError(cudaMalloc(&spts, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cuchebZpoints(degree+1, (cuDoubleComplex*)a, (cuDoubleComplex*)b, spts, 1),__FILE__,__LINE__);
	
		// compute funvals
		cuDoubleComplex *sfvs;
		cuchebCheckError(cudaMalloc(&sfvs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError((*fun)(degree+1, spts, 1, sfvs, 1, USERDATA),__FILE__,__LINE__);
		
		// compute chebcoeffs
		cuDoubleComplex *scfs;
		cuchebCheckError(cudaMalloc(&scfs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cuchebZcoeffs(degree+1, sfvs, 1, scfs, 1),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, scfs, (degree+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		cuchebCheckError(cudaFree(spts),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(scfs),__FILE__,__LINE__);
	}
}

/* user specified tolerance */
ChebPoly::ChebPoly(cuchebCuDoubleComplexFun fun, cuDoubleComplex *A, cuDoubleComplex *B, void *USERDATA, double *tol){

	// set field
	field = CUCHEB_FIELD_DOUBLE_COMPLEX;
	
	// check a and b
	if(cuCreal(*A) == cuCreal(*B) && cuCimag(*A) == cuCimag(*B)){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// check tol
	if(*tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// compute chebpoints
	cuDoubleComplex *zpts;
	cuchebCheckError(cudaMalloc(&zpts, (MAX_DOUBLE_DEG+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cuchebZpoints(MAX_DOUBLE_DEG+1, (cuDoubleComplex*)a, (cuDoubleComplex*)b, zpts, 1),__FILE__,__LINE__);
	
	// compute funvals
	cuDoubleComplex *zfvs;
	cuchebCheckError(cudaMalloc(&zfvs, (MAX_DOUBLE_DEG+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError((*fun)(MAX_DOUBLE_DEG+1, zpts, 1, zfvs, 1, USERDATA),__FILE__,__LINE__);
		
	/* compute chebcoeffs */
	// initialize zcfs
	cuDoubleComplex *zcfs;
	cuchebCheckError(cudaMalloc(&zcfs, (MAX_DOUBLE_DEG+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_DOUBLE_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int max_index;
	int start_index = 0;
	bool converged = false;
	double max_abs, current_abs;
	cuDoubleComplex *max_val, *current_val;
	
	// allocate host pointers
	cuchebCheckError((void*)(max_val = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex))),__FILE__,__LINE__);
	cuchebCheckError((void*)(current_val = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex))),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebZcoeffs(current_degree+1, zfvs, stride, zcfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIzamax(cublasHand, current_degree+1, zcfs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(max_val, &zcfs[max_index-1], sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		max_abs = cuCabs(*max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(current_val, &zcfs[ii], sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			current_abs = cuCabs(*current_val);
			
			// check first coeff
			if(current_abs >= (*tol)*max_abs && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(current_abs >= (*tol)*max_abs && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(current_abs >= (*tol)*max_abs && ii > 1){
				degree = current_degree-ii;
				start_index = ii;
				converged = true;
				break;
			}
			// last coeff
			else if(ii == current_degree-1){
				degree = 0;
				start_index = current_degree;
				converged = true;
				break;
			}
		}
		
		// check current_degree
		if(current_degree > MAX_DOUBLE_DEG){
			printf("\nWarning in %s line: %d\n Function could not be resolved to specified tolerance %e, by a %d degree ChebPoly!\n\n",
				__FILE__,__LINE__,*tol,MAX_DOUBLE_DEG);

			degree = MAX_DOUBLE_DEG;
			start_index = 0;
			converged = true;
		}
	}
	// free host pointers
	free(max_val);
	free(current_val);
	
	// free cublas
	cuchebCheckError(cublasDestroy(cublasHand),__FILE__,__LINE__);
	/* end compute chebcoeffs */
		
	// set coeffs
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &zcfs[start_index], (degree+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(zpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(zfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(zcfs),__FILE__,__LINE__);
}

/* default tolerance */
ChebPoly::ChebPoly(cuchebCuDoubleComplexFun fun, cuDoubleComplex *A, cuDoubleComplex *B, void *USERDATA){

	// set field
	field = CUCHEB_FIELD_DOUBLE_COMPLEX;
	
	// check a and b
	if(cuCreal(*A) == cuCreal(*B) && cuCimag(*A) == cuCimag(*B)){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	cuDoubleComplex *zpts;
	cuchebCheckError(cudaMalloc(&zpts, (MAX_DOUBLE_DEG+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cuchebZpoints(MAX_DOUBLE_DEG+1, (cuDoubleComplex*)a, (cuDoubleComplex*)b, zpts, 1),__FILE__,__LINE__);
	
	// compute funvals
	cuDoubleComplex *zfvs;
	cuchebCheckError(cudaMalloc(&zfvs, (MAX_DOUBLE_DEG+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError((*fun)(MAX_DOUBLE_DEG+1, zpts, 1, zfvs, 1, USERDATA),__FILE__,__LINE__);
		
	/* compute chebcoeffs */
	// initialize zcfs
	cuDoubleComplex *zcfs;
	cuchebCheckError(cudaMalloc(&zcfs, (MAX_DOUBLE_DEG+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_DOUBLE_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int max_index;
	int start_index = 0;
	bool converged = false;
	double max_abs, current_abs;
	cuDoubleComplex *max_val, *current_val;
	
	// allocate host pointers
	cuchebCheckError((void*)(max_val = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex))),__FILE__,__LINE__);
	cuchebCheckError((void*)(current_val = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex))),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebZcoeffs(current_degree+1, zfvs, stride, zcfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIzamax(cublasHand, current_degree+1, zcfs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(max_val, &zcfs[max_index-1], sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		max_abs = cuCabs(*max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(current_val, &zcfs[ii], sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			current_abs = cuCabs(*current_val);
			
			// check first coeff
			if(current_abs >= DBL_EPSILON*max_abs && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(current_abs >= DBL_EPSILON*max_abs && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(current_abs >= DBL_EPSILON*max_abs && ii > 1){
				degree = current_degree-ii;
				start_index = ii;
				converged = true;
				break;
			}
			// last coeff
			else if(ii == current_degree-1){
				degree = 0;
				start_index = current_degree;
				converged = true;
				break;
			}
		}
		
		// check current_degree
		if(current_degree > MAX_DOUBLE_DEG){
			printf("\nWarning in %s line: %d\n Function could not be resolved to machine tolerance %e, by a %d degree ChebPoly!\n\n",
				__FILE__,__LINE__,DBL_EPSILON,MAX_DOUBLE_DEG);

			degree = MAX_DOUBLE_DEG;
			start_index = 0;
			converged = true;
		}
	}
	// free host pointers
	free(max_val);
	free(current_val);
	
	// free cublas
	cuchebCheckError(cublasDestroy(cublasHand),__FILE__,__LINE__);
	/* end compute chebcoeffs */
		
	// set coeffs
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &zcfs[start_index], (degree+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(zpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(zfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(zcfs),__FILE__,__LINE__);
}
