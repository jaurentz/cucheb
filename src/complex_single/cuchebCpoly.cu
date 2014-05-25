#include <cucheb.h>

/* complex single precision constructors */
/* fixed degree */
ChebPoly::ChebPoly(cuchebCuComplexFun fun, cuComplex *A, cuComplex *B, void *USERDATA, int Deg){

	// set field
	field = CUCHEB_FIELD_FLOAT_COMPLEX;

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
	if(cuCrealf(*A) == cuCrealf(*B) && cuCimagf(*A) == cuCimagf(*B)){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// degree 0
	if(degree == 0){
		// compute funvals
		cuComplex *sfvs;
		cuchebCheckError(cudaMalloc(&sfvs, sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError((*fun)(1, (cuComplex*)a, 1, sfvs, 1, USERDATA),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, sfvs, sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);	
	}
	
	// degree > 0
	else{
		// compute chebpoints
		cuComplex *spts;
		cuchebCheckError(cudaMalloc(&spts, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cuchebCpoints(degree+1, (cuComplex*)a, (cuComplex*)b, spts, 1),__FILE__,__LINE__);
	
		// compute funvals
		cuComplex *sfvs;
		cuchebCheckError(cudaMalloc(&sfvs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError((*fun)(degree+1, spts, 1, sfvs, 1, USERDATA),__FILE__,__LINE__);
		
		// compute chebcoeffs
		cuComplex *scfs;
		cuchebCheckError(cudaMalloc(&scfs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cuchebCcoeffs(degree+1, sfvs, 1, scfs, 1),__FILE__,__LINE__);
		
		// set coeffs
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(coeffs, scfs, (degree+1)*sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

		// free device memory
		cuchebCheckError(cudaFree(spts),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(sfvs),__FILE__,__LINE__);
		cuchebCheckError(cudaFree(scfs),__FILE__,__LINE__);
	}
}

/* user specified tolerance */
ChebPoly::ChebPoly(cuchebCuComplexFun fun, cuComplex *A, cuComplex *B, void *USERDATA, float *tol){

	// set field
	field = CUCHEB_FIELD_FLOAT_COMPLEX;
	
	// check a and b
	if(cuCrealf(*A) == cuCrealf(*B) && cuCimagf(*A) == cuCimagf(*B)){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// check tol
	if(*tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// compute chebpoints
	cuComplex *cpts;
	cuchebCheckError(cudaMalloc(&cpts, (MAX_FLOAT_DEG+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cuchebCpoints(MAX_FLOAT_DEG+1, (cuComplex*)a, (cuComplex*)b, cpts, 1),__FILE__,__LINE__);
	
	// compute funvals
	cuComplex *cfvs;
	cuchebCheckError(cudaMalloc(&cfvs, (MAX_FLOAT_DEG+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError((*fun)(MAX_FLOAT_DEG+1, cpts, 1, cfvs, 1, USERDATA),__FILE__,__LINE__);
		
	/* compute chebcoeffs */
	// initialize ccfs
	cuComplex *ccfs;
	cuchebCheckError(cudaMalloc(&ccfs, (MAX_FLOAT_DEG+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_FLOAT_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int max_index;
	int start_index = 0;
	bool converged = false;
	float max_abs, current_abs;
	cuComplex *max_val, *current_val;
	
	// allocate host pointers
	cuchebCheckError((void*)(max_val = (cuComplex*)malloc(sizeof(cuComplex))),__FILE__,__LINE__);
	cuchebCheckError((void*)(current_val = (cuComplex*)malloc(sizeof(cuComplex))),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebCcoeffs(current_degree+1, cfvs, stride, ccfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIcamax(cublasHand, current_degree+1, ccfs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(max_val, &ccfs[max_index-1], sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		max_abs = cuCabsf(*max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(current_val, &ccfs[ii], sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			current_abs = cuCabsf(*current_val);
			
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
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &ccfs[start_index], (degree+1)*sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(cpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(cfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ccfs),__FILE__,__LINE__);
}

/* default tolerance */
ChebPoly::ChebPoly(cuchebCuComplexFun fun, cuComplex *A, cuComplex *B, void *USERDATA){

	// set field
	field = CUCHEB_FIELD_FLOAT_COMPLEX;
	
	// check a and b
	if(cuCrealf(*A) == cuCrealf(*B) && cuCimagf(*A) == cuCimagf(*B)){
		fprintf(stderr,"\nIn %s line: %d, a must not = b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(a, A, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b, B, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	cuComplex *cpts;
	cuchebCheckError(cudaMalloc(&cpts, (MAX_FLOAT_DEG+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cuchebCpoints(MAX_FLOAT_DEG+1, (cuComplex*)a, (cuComplex*)b, cpts, 1),__FILE__,__LINE__);
	
	// compute funvals
	cuComplex *cfvs;
	cuchebCheckError(cudaMalloc(&cfvs, (MAX_FLOAT_DEG+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError((*fun)(MAX_FLOAT_DEG+1, cpts, 1, cfvs, 1, USERDATA),__FILE__,__LINE__);
		
	/* compute chebcoeffs */
	// initialize ccfs
	cuComplex *ccfs;
	cuchebCheckError(cudaMalloc(&ccfs, (MAX_FLOAT_DEG+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	
	// initialize compute variables
	int stride = pow(2,MAX_FLOAT_DEG_EXP-3); 
	int current_degree = pow(2,3);
	int max_index;
	int start_index = 0;
	bool converged = false;
	float max_abs, current_abs;
	cuComplex *max_val, *current_val;
	
	// allocate host pointers
	cuchebCheckError((void*)(max_val = (cuComplex*)malloc(sizeof(cuComplex))),__FILE__,__LINE__);
	cuchebCheckError((void*)(current_val = (cuComplex*)malloc(sizeof(cuComplex))),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublasHand;
	cuchebCheckError(cublasCreate(&cublasHand),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublasHand, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// compute coeffs adaptively until convergence
	while(converged != true){
		// compute cheb interpolant of current_degree 
		cuchebCheckError(cuchebCcoeffs(current_degree+1, cfvs, stride, ccfs, 1),__FILE__,__LINE__);

		// get max_index
		cuchebCheckError(cublasIcamax(cublasHand, current_degree+1, ccfs, 1, &max_index),__FILE__,__LINE__);
		
		// set maximum modulus of coefficient
		cuchebCheckError(cudaMemcpy(max_val, &ccfs[max_index-1], sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		max_abs = cuCabsf(*max_val);

		// check for convergence
		for(int ii=0;ii<current_degree;ii++){
			// get current coefficient
			cuchebCheckError(cudaMemcpy(current_val, &ccfs[ii], sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
			current_abs = cuCabsf(*current_val);
			
			// check first coeff
			if(current_abs >= FLT_EPSILON*max_abs && ii == 0){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check second coeff
			else if(current_abs >= FLT_EPSILON*max_abs && ii == 1){
				stride = stride/2;
				current_degree = current_degree*2;
				converged = false;
				break;
			}
			// check middle coeffs
			else if(current_abs >= FLT_EPSILON*max_abs && ii > 1){
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
	cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(coeffs, &ccfs[start_index], (degree+1)*sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free device memory
	cuchebCheckError(cudaFree(cpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(cfvs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ccfs),__FILE__,__LINE__);
}
