#include <cucheb.h>

/* double precision constructors */
/* fixed degree */
ChebPoly::ChebPoly(cuchebDoubleFun fun, void *USERDATA, double *A, double *B, int Deg){

	// set field
	field = CUCHEB_FIELD_DOUBLE;

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
	if(*A == *B){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// degree 0
	if(degree == 0){
		// compute funvals
		double left, right;
		cuchebCheckError((*fun)(1, A, 1, &left, 1, USERDATA),__FILE__,__LINE__);
		cuchebCheckError((*fun)(1, B, 1, &right, 1, USERDATA),__FILE__,__LINE__);
		
		// set coeffs
		left = (right+left)/2.0;
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, &left, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);

	}
	
	// degree > 0
	else{
		// compute chebpoints
		double *pts;
		cuchebCheckError((void*)(pts = (double*)malloc((degree+1)*sizeof(double))),__FILE__,__LINE__);
		cuchebCheckError(cuchebDpoints(degree+1, A, B, pts, 1),__FILE__,__LINE__);
	
		// compute funvals
		double *fvs;
		cuchebCheckError((void*)(fvs = (double*)malloc((degree+1)*sizeof(double))),__FILE__,__LINE__);
		cuchebCheckError((*fun)(degree+1, pts, 1, fvs, 1, USERDATA),__FILE__,__LINE__);
		double *dfvs;
		cuchebCheckError(cudaMalloc(&dfvs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(dfvs, fvs, (degree+1)*sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
		
		// compute chebcoeffs
		double *dcfs;
		cuchebCheckError(cudaMalloc(&dcfs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cuchebDcoeffs(degree+1, dfvs, 1, dcfs, 1),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, dcfs, (degree+1)*sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		free(pts);
		free(fvs);
		cuchebCheckError(cudaFree(dfvs),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(dcfs),__FILE__,__LINE__);
	}
}

/* user specified tolerance */
ChebPoly::ChebPoly(cuchebDoubleFun fun, void *USERDATA, double *A, double *B, double *tol){

	// set field
	field = CUCHEB_FIELD_DOUBLE;
	
	// check a and b
	if(*A == *B){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// check tol
	if(*tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// compute chebpoints
	double *pts;
	cuchebCheckError((void*)(pts = (double*)malloc((MAX_DOUBLE_DEG+1)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cuchebDpoints(MAX_DOUBLE_DEG+1, A, B, pts, 1),__FILE__,__LINE__);
	
	// compute funvals
	double *fvs;
	cuchebCheckError((void*)(fvs = (double*)malloc((MAX_DOUBLE_DEG+1)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError((*fun)(MAX_DOUBLE_DEG+1, pts, 1, fvs, 1, USERDATA),__FILE__,__LINE__);
	double *dfvs;
	cuchebCheckError(cudaMalloc(&dfvs, (MAX_DOUBLE_DEG+1)*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dfvs, fvs, (MAX_DOUBLE_DEG+1)*sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// get max_index
	int max_index;
	cuchebCheckError(cublasIdamax(cublasHand, MAX_DOUBLE_DEG+1, dfvs, 1, &max_index),__FILE__,__LINE__);
		
	// set maximum modulus of coefficient
	double max_val;
	cuchebCheckError(cudaMemcpy(&max_val, &dfvs[max_index-1], sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	max_val = abs(max_val);
		
	/* compute chebcoeffs */
	// initialize dcfs
	double *dcfs;
	cuchebCheckError(cudaMalloc(&dcfs, (MAX_DOUBLE_DEG+1)*sizeof(double)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_DOUBLE_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int start_index = 0;
	bool converged = false;
	double current_val;
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebDcoeffs(current_degree+1, dfvs, stride, dcfs, 1),__FILE__,__LINE__);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(&current_val, &dcfs[ii], sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			current_val = abs(current_val);
			
			// check first coeff
			if(current_val >= (*tol)*max_val && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(current_val >= (*tol)*max_val && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(current_val >= (*tol)*max_val && ii > 1){
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
	free(pts);
	free(fvs);
	
	// free cublas
	cuchebCheckError(cublasDestroy(cublasHand),__FILE__,__LINE__);
	/* end compute chebcoeffs */
		
	// set coeffs
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &dcfs[start_index], (degree+1)*sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(dfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dcfs),__FILE__,__LINE__);
}

/* default tolerance */
ChebPoly::ChebPoly(cuchebDoubleFun fun, void *USERDATA, double *A, double *B){

	// set field
	field = CUCHEB_FIELD_DOUBLE;
	
	// check a and b
	if(*A == *B){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	double *pts;
	cuchebCheckError((void*)(pts = (double*)malloc((MAX_DOUBLE_DEG+1)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cuchebDpoints(MAX_DOUBLE_DEG+1, A, B, pts, 1),__FILE__,__LINE__);
	
	// compute funvals
	double *fvs;
	cuchebCheckError((void*)(fvs = (double*)malloc((MAX_DOUBLE_DEG+1)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError((*fun)(MAX_DOUBLE_DEG+1, pts, 1, fvs, 1, USERDATA),__FILE__,__LINE__);
	double *dfvs;
	cuchebCheckError(cudaMalloc(&dfvs, (MAX_DOUBLE_DEG+1)*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dfvs, fvs, (MAX_DOUBLE_DEG+1)*sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);

	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// get max_index
	int max_index;
	cuchebCheckError(cublasIdamax(cublasHand, MAX_DOUBLE_DEG+1, dfvs, 1, &max_index),__FILE__,__LINE__);
		
	// set maximum modulus of coefficient
	double max_val;
	cuchebCheckError(cudaMemcpy(&max_val, &dfvs[max_index-1], sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	max_val = abs(max_val);
	
	/* compute chebcoeffs */
	// initialize dcfs
	double *dcfs;
	cuchebCheckError(cudaMalloc(&dcfs, (MAX_DOUBLE_DEG+1)*sizeof(double)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_DOUBLE_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int start_index = 0;
	bool converged = false;
	double current_val;
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebDcoeffs(current_degree+1, dfvs, stride, dcfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIdamax(cublasHand, current_degree+1, dfvs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(&max_val, &dfvs[max_index-1], sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		max_val = abs(max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(&current_val, &dcfs[ii], sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			current_val = abs(current_val);
			
			// check first coeff
			if(current_val >= DBL_EPSILON*max_val && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(current_val >= DBL_EPSILON*max_val && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(current_val >= DBL_EPSILON*max_val && ii > 1){
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
	free(pts);
	free(fvs);
	
	// free cublas
	cuchebCheckError(cublasDestroy(cublasHand),__FILE__,__LINE__);
	/* end compute chebcoeffs */
		
	// set coeffs
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &dcfs[start_index], (degree+1)*sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(dfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dcfs),__FILE__,__LINE__);
}
