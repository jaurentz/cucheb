#include <cucheb.h>

cuchebStatus_t cuchebDeigs(cuchebLanczosHandle* LH, cuchebOpMult OPMULT, void* USERDATA, double *eigvecs){
	// check n
	int n = LH->n;
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check numeigs
	int numeigs = LH->numeigs;
	if(numeigs < 1 || numeigs > n){
		fprintf(stderr,"\nIn %s line: %d, numeigs must be > 0 and <= n.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check runlength
	int runlength = LH->runlength;
	if(runlength < 1 || runlength >= n){
		fprintf(stderr,"\nIn %s line: %d, runlength must be > 0 and < n.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check restarts
	int restarts = LH->restarts;
	if(restarts < 0){
		fprintf(stderr,"\nIn %s line: %d, restarts must be => 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check tol
	double tol = LH->tol;
	if(tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(tol < DBL_EPSILON){
		fprintf(stderr,"\nIn %s line: %d, tol is below machine precision. Algorithm may not converge.\n",__FILE__,__LINE__);
	}

	// check numconv
	int numconv = LH->numconv;
	if(numconv < 0){numconv = 0;}
	if(numconv >= numeigs){return CUCHEB_STATUS_SUCCESS;}
	
	// allocate memory for Lanzcos
	double *vecs, *diags, *sdiags, *ritzvecs;
	cuchebCheckError(cudaMalloc(&vecs,(runlength+1)*n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&diags,runlength*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&sdiags,runlength*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ritzvecs,(numeigs+1)*n*sizeof(double)),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// check starting vector
	double temp;
	curandGenerator_t curand_gen;
	cuchebCheckError(cublasDnrm2(cublas_handle,n,eigvecs,1,&temp),__FILE__,__LINE__);
	if(temp >= 0.0){
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,eigvecs,1),__FILE__,__LINE__);
	}
	else{
		// initialize curand
		cuchebCheckError(curandCreateGenerator(&curand_gen, CURAND_RNG_PSEUDO_DEFAULT),__FILE__,__LINE__);
		cuchebCheckError(curandSetPseudoRandomGeneratorSeed(curand_gen,time(NULL)),__FILE__,__LINE__);
	
		// random starting vector
		cuchebCheckError(curandGenerateNormalDouble(curand_gen,vecs,n,0.0,1.0),__FILE__,__LINE__);
		cuchebCheckError(cublasDnrm2(cublas_handle,n,vecs,1,&temp),__FILE__,__LINE__);
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,vecs,1),__FILE__,__LINE__);	
		
		// shutdown curand
		cuchebCheckError(curandDestroyGenerator(curand_gen),__FILE__,__LINE__);
	}
	cuchebCheckError(cudaMemcpy(vecs,eigvecs,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	// lanczos
	int nummatvecs = 0;
	for(int ii=0;ii<restarts+1;ii++){
	
		// do Lanzcos run
		cuchebCheckError(cuchebDlanczos(n,OPMULT,USERDATA,numconv,runlength,vecs,diags,sdiags),__FILE__,__LINE__);
		
		// update nummatvecs
		nummatvecs += runlength-numconv;

		// restart
		cuchebCheckError(cuchebDrestart(n,runlength,numeigs,&numconv,vecs,diags,sdiags,ritzvecs,tol),__FILE__,__LINE__);

		// check convergence
		if(numconv == numeigs){
			LH->numconv = numconv;
			LH->numrestarts = ii;
			LH->nummatvecs = nummatvecs;
			break;
		}
		
		// check iterations
		if(ii == restarts){
			LH->numconv = numconv;
			LH->numrestarts = ii;
			LH->nummatvecs = nummatvecs;
		}
	}	
			
	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);
	
	// copy ritzvecs into eigvecs
	cuchebCheckError(cudaMemcpy(eigvecs,ritzvecs,numeigs*n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

	// free memory
	cuchebCheckError(cudaFree(vecs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(diags),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(sdiags),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ritzvecs),__FILE__,__LINE__);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}

cuchebStatus_t cuchebDeigs(cuchebLanczosHandle* LH, ChebOp* CO, double *eigvecs){

	// check n
	int n = LH->n;
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check numeigs
	int numeigs = LH->numeigs;
	if(numeigs < 1 || numeigs > n){
		fprintf(stderr,"\nIn %s line: %d, numeigs must be > 0 and <= n.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check runlength
	int runlength = LH->runlength;
	if(runlength < 1 || runlength >= n){
		fprintf(stderr,"\nIn %s line: %d, runlength must be > 0 and < n.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check restarts
	int restarts = LH->restarts;
	if(restarts < 0){
		fprintf(stderr,"\nIn %s line: %d, restarts must be => 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check tol
	double tol = LH->tol;
	if(tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(tol < DBL_EPSILON){
		fprintf(stderr,"\nIn %s line: %d, tol is below machine precision. Algorithm may not converge.\n",__FILE__,__LINE__);
	}

	// check numconv
	int numconv = LH->numconv;
	if(numconv < 0){numconv = 0;}
	if(numconv >= numeigs){return CUCHEB_STATUS_SUCCESS;}

	// allocate memory for Lanzcos
	double *vecs, *diags, *sdiags, *ritzvecs;
	cuchebCheckError(cudaMalloc(&vecs,(runlength+1)*n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&diags,runlength*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&sdiags,runlength*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ritzvecs,(numeigs+1)*n*sizeof(double)),__FILE__,__LINE__);
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);
	
	// check starting vector
	double temp;
	curandGenerator_t curand_gen;
	cuchebCheckError(cublasDnrm2(cublas_handle,n,eigvecs,1,&temp),__FILE__,__LINE__);
	if(temp >= 0.0){
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,eigvecs,1),__FILE__,__LINE__);
	}
	else{
		// initialize curand
		cuchebCheckError(curandCreateGenerator(&curand_gen, CURAND_RNG_PSEUDO_DEFAULT),__FILE__,__LINE__);
		cuchebCheckError(curandSetPseudoRandomGeneratorSeed(curand_gen,time(NULL)),__FILE__,__LINE__);
	
		// random starting vector
		cuchebCheckError(curandGenerateNormalDouble(curand_gen,vecs,n,0.0,1.0),__FILE__,__LINE__);
		cuchebCheckError(cublasDnrm2(cublas_handle,n,vecs,1,&temp),__FILE__,__LINE__);
		temp = 1.0/temp;
		cuchebCheckError(cublasDscal(cublas_handle,n,&temp,vecs,1),__FILE__,__LINE__);	
		
		// shutdown curand
		cuchebCheckError(curandDestroyGenerator(curand_gen),__FILE__,__LINE__);
	}
	cuchebCheckError(cudaMemcpy(vecs,eigvecs,n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	
	// lanczos
	int nummatvecs = 0;
	for(int ii=0;ii<restarts+1;ii++){
	
		// do Lanzcos run
		cuchebCheckError(cuchebDlanczos(CO,numconv,runlength,vecs,diags,sdiags),__FILE__,__LINE__);
		
		// update nummatvecs
		nummatvecs += (runlength-numconv)*(CO->getChebpoly()->getDegree());

		// restart
		cuchebCheckError(cuchebDrestart(n,runlength,numeigs,&numconv,vecs,diags,sdiags,ritzvecs,tol),__FILE__,__LINE__);

		// check convergence
		if(numconv == numeigs){
			LH->numconv = numconv;
			LH->numrestarts = ii;
			LH->nummatvecs = nummatvecs;
			break;
		}
		
		// check iterations
		if(ii == restarts){
			LH->numconv = numconv;
			LH->numrestarts = ii;
			LH->nummatvecs = nummatvecs;
		}
	}	
			
	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);
	
	// copy ritzvecs into eigvecs
	cuchebCheckError(cudaMemcpy(eigvecs,ritzvecs,numeigs*n*sizeof(double),cudaMemcpyDeviceToDevice),__FILE__,__LINE__);

	// free memory
	cuchebCheckError(cudaFree(vecs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(diags),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(sdiags),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ritzvecs),__FILE__,__LINE__);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}
