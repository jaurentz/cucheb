#include <cucheb.h>

/* single precision constructors */
/* fixed degree */
ChebPoly::ChebPoly(cuchebFloatFun fun, float *A, float *B, void *USERDATA, int Deg){

	// set field
	field = CUCHEB_FIELD_FLOAT;

	// check degree
	if(Deg < 0){
		fprintf(stderr,"\nIn %s line: %d, degree must be >= 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(Deg > MAX_FLOAT_DEG){
		fprintf(stderr,"\nIn %s line: %d, degree must be <= %d.\n",__FILE__,__LINE__,MAX_FLOAT_DEG);
		cuchebExit(-1);
	}
	
	// set degree
	degree = Deg;
	
	// check a and b
	if(*A == *B){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// degree 0
	if(degree == 0){
		// compute funvals
		float *sfvs;
		cuchebCheckError(cudaMalloc(&sfvs, sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError((*fun)(1, (float*)a, 1, sfvs, 1,USERDATA),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, sfvs, sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);	
	}
	
	// degree > 0
	else{
		// compute chebpoints
		float *spts;
		cuchebCheckError(cudaMalloc(&spts, (degree+1)*sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cuchebSpoints(degree+1, (float*)a, (float*)b, spts, 1),__FILE__,__LINE__);
	
		// compute funvals
		float *sfvs;
		cuchebCheckError(cudaMalloc(&sfvs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError((*fun)(degree+1, spts, 1, sfvs, 1,USERDATA),__FILE__,__LINE__);
		
		// compute chebcoeffs
		float *scfs;
		cuchebCheckError(cudaMalloc(&scfs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cuchebScoeffs(degree+1, sfvs, 1, scfs, 1),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, scfs, (degree+1)*sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		cuchebCheckError(cudaFree(spts),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(scfs),__FILE__,__LINE__);
	}
}

/* user specified tolerance */
ChebPoly::ChebPoly(cuchebFloatFun fun, float *A, float *B, void *USERDATA, float *tol){

	// set field
	field = CUCHEB_FIELD_FLOAT;
	
	// check a and b
	if(*A == *B){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// check tol
	if(*tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// compute chebpoints
	float *spts;
	cuchebCheckError(cudaMalloc(&spts, (MAX_FLOAT_DEG+1)*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cuchebSpoints((MAX_FLOAT_DEG+1), (float*)a, (float*)b, spts, 1),__FILE__,__LINE__);
	
	// compute funvals
	float *sfvs;
	cuchebCheckError(cudaMalloc(&sfvs, (MAX_FLOAT_DEG+1)*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError((*fun)((MAX_FLOAT_DEG+1), spts, 1, sfvs, 1,USERDATA),__FILE__,__LINE__);
		
	/* compute chebcoeffs */
	// initialize scfs
	float *scfs;
	cuchebCheckError(cudaMalloc(&scfs, (MAX_FLOAT_DEG+1)*sizeof(float)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_FLOAT_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int max_index;
	int start_index = 0;
	bool converged = false;
	float *max_val, *current_val;
	
	// allocate host pointers
	cuchebCheckError((void*)(max_val = (float*)malloc(sizeof(float))),__FILE__,__LINE__);
	cuchebCheckError((void*)(current_val = (float*)malloc(sizeof(float))),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebScoeffs(current_degree+1, sfvs, stride, scfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIsamax(cublasHand, current_degree+1, scfs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(max_val, &scfs[max_index-1], sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		*max_val = abs(*max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(current_val, &scfs[ii], sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			*current_val = abs(*current_val);
			
			// check first coeff
			if(*current_val >= (*tol)*(*max_val) && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(*current_val >= (*tol)*(*max_val) && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(*current_val >= (*tol)*(*max_val) && ii > 1){
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
		if(current_degree > MAX_FLOAT_DEG){
			printf("\nWarning in %s line: %d\n Function could not be resolved to specified tolerance %e, by a %d degree ChebPoly!\n\n",
				__FILE__,__LINE__,*tol,MAX_FLOAT_DEG);

			degree = MAX_FLOAT_DEG;
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
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &scfs[start_index], (degree+1)*sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(spts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(scfs),__FILE__,__LINE__);
}

/* default tolerance */
ChebPoly::ChebPoly(cuchebFloatFun fun, float *A, float *B, void *USERDATA){

	// set field
	field = CUCHEB_FIELD_FLOAT;
	
	// check a and b
	if(*A == *B){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	float *spts;
	cuchebCheckError(cudaMalloc(&spts, (MAX_FLOAT_DEG+1)*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cuchebSpoints((MAX_FLOAT_DEG+1), (float*)a, (float*)b, spts, 1),__FILE__,__LINE__);
	
	// compute funvals
	float *sfvs;
	cuchebCheckError(cudaMalloc(&sfvs, (MAX_FLOAT_DEG+1)*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError((*fun)((MAX_FLOAT_DEG+1), spts, 1, sfvs, 1,USERDATA),__FILE__,__LINE__);
		
	/* compute chebcoeffs */
	// initialize scfs
	float *scfs;
	cuchebCheckError(cudaMalloc(&scfs, (MAX_FLOAT_DEG+1)*sizeof(float)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_FLOAT_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int max_index;
	int start_index = 0;
	bool converged = false;
	float *max_val, *current_val;
	
	// allocate host pointers
	cuchebCheckError((void*)(max_val = (float*)malloc(sizeof(float))),__FILE__,__LINE__);
	cuchebCheckError((void*)(current_val = (float*)malloc(sizeof(float))),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebScoeffs(current_degree+1, sfvs, stride, scfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIsamax(cublasHand, current_degree+1, scfs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(max_val, &scfs[max_index-1], sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		*max_val = abs(*max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(current_val, &scfs[ii], sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			*current_val = abs(*current_val);
			
			// check first coeff
			if(*current_val >= FLT_EPSILON*(*max_val) && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(*current_val >= FLT_EPSILON*(*max_val) && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(*current_val >= FLT_EPSILON*(*max_val) && ii > 1){
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
		if(current_degree > MAX_FLOAT_DEG){
			printf("\nWarning in %s line: %d\n Function could not be resolved to machine tolerance %e, by a %d degree ChebPoly!\n\n",
				__FILE__,__LINE__,FLT_EPSILON,MAX_FLOAT_DEG);

			degree = MAX_FLOAT_DEG;
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
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &scfs[start_index], (degree+1)*sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(spts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(scfs),__FILE__,__LINE__);
}

