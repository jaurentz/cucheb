#include <cucheb.h>

__global__ void	DcudaEliminateKernel(double *a, double *b, double *d_diag, double *d_sdiag, double *d_V);
void DcudaEliminate(double *a, double *b, double *d_diag, double *d_sdiag, double *d_V);
__global__ void	DcudaRotMultKernel(int N, double *d_Q, double *d_V);
void DcudaRotMult(int N, double *d_Q, double *d_V);
__global__ void	DcudaBulgeKernel(double *d_sdiag, double *d_V); 
void DcudaBulge(double *d_sdiag, double *d_V);
__global__ void	DcudaTwiddleKernel(double *d_sdiag, double *d_h, double *d_Z);
void DcudaTwiddle(double *d_sdiag, double *d_h, double *d_Z);

/* restart routine for Lanczos */
cuchebStatus_t cuchebDrestart(int n,int runlength,int neigs,int *nconv,double *vecs,double *diags,double *sdiags,double *ritzvecs,double tol){

	// check n
	if(n < 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check neigs
	if(neigs < 1){
		fprintf(stderr,"\nIn %s line: %d, neigs must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check runlength
	if(runlength <= neigs){
		fprintf(stderr,"\nIn %s line: %d, runlength must be > neigs.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check nconv
	if((*nconv) < 0){
		fprintf(stderr,"\nIn %s line: %d, nconv must be => 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(neigs <= (*nconv)){
		fprintf(stderr,"\nIn %s line: %d, nconv must be < neigs.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check tol
	if(tol <= 0){
		fprintf(stderr,"\nIn %s line: %d, tol must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// initialize cublas
	cublasHandle_t cublas_handle;
	cuchebCheckError(cublasCreate(&cublas_handle),__FILE__,__LINE__);
	cuchebCheckError(cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST),__FILE__,__LINE__);

	// initialize gpu eigenvectors
	int len = runlength-*nconv;
	double temp;
	double *eigvecs;
	cuchebCheckError(cudaMalloc(&eigvecs,len*len*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(len*len,eigvecs,1,0.0),__FILE__,__LINE__);
	temp = 1.0;
	for(int ii=0;ii<len;ii++){
		cuchebCheckError(cublasSetVector(1,sizeof(double),&temp,1,&eigvecs[ii+ii*len],1),__FILE__,__LINE__);
	}

	// initialize cpu diags
	double *h_diags;
	cuchebCheckError((void*)(h_diags = (double*)malloc(len*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(h_diags,diags,len*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

	// initialize cpu sdiags
	double *h_sdiags;
	cuchebCheckError((void*)(h_sdiags = (double*)malloc((len-1)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(h_sdiags,sdiags,(len-1)*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

	// initialize cpu eigenvectors
	double *h_eigvecs;
	cuchebCheckError((void*)(h_eigvecs = (double*)malloc(len*len*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(h_eigvecs,eigvecs,len*len*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);

	// call lapacke
	lapack_int lpint;
	lpint = LAPACKE_dsteqr(LAPACK_COL_MAJOR,'I',len,h_diags,h_sdiags,h_eigvecs,len);
	if(lpint |= 0){
		fprintf(stderr,"\nLapacke error: %d occured in %s at line: %d\n\n",lpint,__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// copy results back to gpu
	cuchebCheckError(cudaMemcpy(diags,h_diags,len*sizeof(double),cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(eigvecs,h_eigvecs,len*len*sizeof(double),cudaMemcpyHostToDevice),__FILE__,__LINE__);

	// swap eigenvalues
	for(int ii=0;ii < (len/2);ii++){
		cuchebCheckError(cublasDswap(cublas_handle,1,&diags[ii],1,&diags[len-ii-1],1),__FILE__,__LINE__);
	}
	
//for(int jj=0;jj<runlength;jj++){
//	cudaMemcpy(&temp, &diags[jj], sizeof(double), cudaMemcpyDeviceToHost);
//	printf("diags[%d] = %e\n",jj,temp);
//}
//for(int jj=0;jj<runlength;jj++){
//	cudaMemcpy(&temp, &sdiags[jj], sizeof(double), cudaMemcpyDeviceToHost);
//	printf("sdiags[%d] = %e\n",jj,temp);
//}
//for(int jj=0;jj<len;jj++){
//	cudaMemcpy(&temp, &eigvecs[(jj+1)*(len)-1], sizeof(double), cudaMemcpyDeviceToHost);
//	printf("eigvecs[%d] = %e\n",(jj+1)*(len)-1,temp);
//}	
	// swap eigenvectors
	for(int ii=0;ii < (len/2);ii++){
		cuchebCheckError(cublasDswap(cublas_handle,len,&eigvecs[ii*len],1,&eigvecs[(len-ii-1)*len],1),__FILE__,__LINE__);
	}
	
//for(int jj=0;jj<8;jj++){
//	cudaMemcpy(&temp, &vecs[jj], sizeof(double), cudaMemcpyDeviceToHost);
//	printf("vecs[%d] = %e\n",jj,temp);
//}
                           
	// compute ritz vecs
    temp = 1.0;
    double temp2 = 0.0;
	cuchebCheckError(cublasDgemm(cublas_handle,CUBLAS_OP_N,CUBLAS_OP_N,n,(neigs-(*nconv)+1),len,&temp,&vecs[(*nconv)*n],n,
		eigvecs,len,&temp2,ritzvecs,n),__FILE__,__LINE__);

	// copy neigs-nconv+1 ritzvecs into vecs
	cuchebCheckError(cudaMemcpy(&vecs[(*nconv)*n],ritzvecs,(neigs-(*nconv)+1)*n*sizeof(double),
		cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		
//for(int jj=0;jj<8;jj++){
//	cudaMemcpy(&temp, &vecs[jj], sizeof(double), cudaMemcpyDeviceToHost);
//	printf("vecs[%d] = %e\n",jj,temp);
//}

	// check convergence
	double temp3;
	int neweigs = 0;
	cuchebCheckError(cublasDnrm2(cublas_handle,1,&sdiags[runlength-1],1,&temp),__FILE__,__LINE__);
//printf("\ntemp = %e\n",temp);
	cuchebCheckError(cublasDnrm2(cublas_handle,1,diags,1,&temp3),__FILE__,__LINE__);
//printf("temp3 = %e\n",temp3);
//printf("tol = %e\n",tol);
	for(int ii=0;ii<(neigs-(*nconv));ii++){
		cuchebCheckError(cublasDnrm2(cublas_handle,1,&eigvecs[(ii+1)*(len)-1],1,&temp2),__FILE__,__LINE__);
//printf("temp2 = %e\n",temp2);
		// lock converged
		if(temp*temp2 < tol*temp3){
			temp2 = 0.0;
			cuchebCheckError(cublasSetVector(1,sizeof(double),&temp2,1,&eigvecs[(ii+1)*(len)-1],1),__FILE__,__LINE__);
			neweigs += 1;
		}
	}
//printf("nconv = %d\n",*nconv);
//printf("neweigs = %d\n",neweigs);		
	// update nconv
	*nconv += neweigs;

	// exit if nconv = neigs
	if((*nconv) == neigs){

		// shutdown cublas
		cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);

		// free memory
		cuchebCheckError(cudaFree(eigvecs),__FILE__,__LINE__);

		free(h_eigvecs);
		free(h_diags);
		free(h_sdiags);

		// return success
		return CUCHEB_STATUS_SUCCESS;
	}

//cudaMemcpy(&temp, &eigvecs[(*nconv+1)*(len)-1], sizeof(double), cudaMemcpyDeviceToHost);
//printf("temp = %e\n",temp);
//cudaMemcpy(&temp, &eigvecs[(*nconv+2)*(len)-1], sizeof(double), cudaMemcpyDeviceToHost);
//printf("temp = %e\n",temp);
//cudaMemcpy(&temp, &diags[(*nconv)], sizeof(double), cudaMemcpyDeviceToHost);
//printf("temp = %e\n",temp);
//cudaMemcpy(&temp, &diags[(*nconv+1)], sizeof(double), cudaMemcpyDeviceToHost);
//printf("temp = %e\n",temp);
//cudaMemcpy(&temp, &sdiags[(*nconv)], sizeof(double), cudaMemcpyDeviceToHost);
//printf("temp = %e\n",temp);

	// update restart vectors
	for(int ii=*nconv;ii<neigs;ii++){
		// create zero in d_Z
		DcudaEliminate(&eigvecs[(ii+1)*(len)-1],&eigvecs[(ii+2)*(len)-1],&diags[ii],&sdiags[ii],ritzvecs);

		// update vecs
		DcudaRotMult(n,&vecs[ii*n],ritzvecs);

		// chase bulge
		for(int jj=0;jj<ii-*nconv;jj++){
			// build bulge
			DcudaBulge(&sdiags[ii-jj-1],ritzvecs);

			// remove bulge
			DcudaEliminate(&ritzvecs[2],&sdiags[ii-jj],&diags[ii-jj-1],&sdiags[ii-jj-1],ritzvecs);

			// update vecs
			DcudaRotMult(n,&vecs[(ii-jj-1)*n],ritzvecs);
		}
	}

	// update last entry of s_diag
	DcudaTwiddle(&sdiags[neigs],&sdiags[runlength-1],&eigvecs[(neigs+1)*len-1]);

	// shutdown cublas
	cuchebCheckError(cublasDestroy(cublas_handle),__FILE__,__LINE__);

	// free memory
	cuchebCheckError(cudaFree(eigvecs),__FILE__,__LINE__);

	free(h_eigvecs);
	free(h_diags);
	free(h_sdiags);

	// return success
	return CUCHEB_STATUS_SUCCESS;
}

__global__ void	DcudaEliminateKernel(double *a, double *b, double *d_diag, double *d_sdiag, double *d_V){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	double temp1, temp2, temp3;
	double cc, cs, ss;
	
	if(ii < 1){
		// build rotator generators
		if(abs(*a) == 0){
			d_V[0] = 1.0;
			d_V[1] = 0.0;

			// update d_Z
			(*a) = 0.0f;
		}
		else if(abs(*a) > abs(*b)){
			temp2 = *b/(*a);
			temp1 = abs(*a)*sqrt(1.0 + temp2*temp2);
			d_V[0] = (*b)/temp1;
			d_V[1] = (*a)/temp1;

			// update d_Z
			(*a) = 0.0f;
			(*b) = temp1;
		}
		else{
			temp2 = *a/(*b);
			temp1 = abs(*b)*sqrt(1.0 + temp2*temp2);
			d_V[0] = (*b)/temp1;
			d_V[1] = (*a)/temp1;

			// update d_Z
			(*a) = 0.0;
			(*b) = temp1;
		}

		// ensure orthogonality
		temp1 = sqrt(d_V[0]*d_V[0] + d_V[1]*d_V[1]);
		d_V[0] = d_V[0]/temp1;
		d_V[1] = d_V[1]/temp1;

		// update d_diag and d_sdiag
		cc = d_V[0]*d_V[0];
		cs = d_V[1]*d_V[0];
		ss = d_V[1]*d_V[1];
		temp1 = cc*d_diag[0] - 2.0*cs*d_sdiag[0] + ss*d_diag[1];
		temp2 = cs*d_diag[0] + (cc-ss)*d_sdiag[0] - cs*d_diag[1];
		temp3 = ss*d_diag[0] + 2.0*cs*d_sdiag[0] + cc*d_diag[1];
		d_diag[0] = temp1;
		d_diag[1] = temp3;
		d_sdiag[0] = temp2;	
	}
}
void DcudaEliminate(double *a, double *b, double *d_diag, double *d_sdiag, double *d_V){
	// call kernel
	DcudaEliminateKernel<<<1,1>>>(a,b,d_diag,d_sdiag,d_V);
} 

__global__ void	DcudaRotMultKernel(int N, double *d_Q, double *d_V){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;
	double temp1;
	
	if(ii < N){
		// update d_Q
		temp1 = d_V[0]*d_Q[ii] - d_V[1]*d_Q[N+ii];
		d_Q[N+ii] = d_V[1]*d_Q[ii] + d_V[0]*d_Q[N+ii];
		d_Q[ii] = temp1;		
	}
}
void DcudaRotMult(int N, double *d_Q, double *d_V){
	// compute variables
	int BLOCK_SIZE, GRID_SIZE;

	// set BLOCK_SIZE and compute GRID_SIZE
	BLOCK_SIZE = 512;
	if(N%BLOCK_SIZE == 0){GRID_SIZE = N/512;}
	else{GRID_SIZE = N/512 + 1;}

	// call kernel
	DcudaRotMultKernel<<<GRID_SIZE,BLOCK_SIZE>>>(N,d_Q,d_V);
}
__global__ void	DcudaBulgeKernel(double *d_sdiag, double *d_V){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < 1){
		// build bulge
		d_V[2] = d_V[1]*d_sdiag[0];

		// update d_sdiag
		d_sdiag[0] = d_V[0]*d_sdiag[0];		
	}
}  
void DcudaBulge(double *d_sdiag, double *d_V){
	// call kernel
	DcudaBulgeKernel<<<1,1>>>(d_sdiag,d_V);
}
__global__ void	DcudaTwiddleKernel(double *d_sdiag, double *d_h, double *d_Z){
	int tix = threadIdx.x, bix = blockIdx.x, bdx = blockDim.x;
	int ii = bix*bdx+tix;

	if(ii < 1){
		// adjust twiddle
		d_sdiag[0] = d_h[0]*d_Z[0];	
	}
} 
void DcudaTwiddle(double *d_sdiag, double *d_h, double *d_Z){
	// call kernel
	DcudaTwiddleKernel<<<1,1>>>(d_sdiag,d_h,d_Z);
}

