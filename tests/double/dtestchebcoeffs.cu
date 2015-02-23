#include <cucheb.h>

/* helpers for double precision constructors */
cuchebStatus_t dfuncaller(int n, const double *in, int incin, double *out, int incout, void* userdata){

	double *tau = (double*)userdata;
	
	for(int ii=0;ii<n;ii++){
		out[ii*incout] = sin((*tau)*in[ii*incin]);
	}
	
	return CUCHEB_STATUS_SUCCESS;

}

// driver
int main(){

	// compute variables
	/* double */
	int dn = pow(2,5)+1;
	int sizeD = sizeof(double);
	double *dpts, *dfvs, *dcfs;
	double da, db; 
	double *ddcfs, *ddfvs;

        // set device
        cuchebCheckError(cudaSetDevice(7),__FILE__,__LINE__); 
	
	// allocate host memory
	/* double */
	cuchebCheckError((void*)(dpts = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(dcfs = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(dfvs = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	
	// allocate device memory
	/* double */
	cuchebCheckError(cudaMalloc(&ddcfs, dn*sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ddfvs, dn*sizeD),__FILE__,__LINE__);
	
	// set host pointers
	/* double */
	da = -1.0;
	db = 1.0;
	
	// compute chebpoints
	/* double */
	cuchebCheckError(cuchebDpoints(dn,&da,&db,dpts,1),__FILE__,__LINE__);
	
	// compute funvals
	/* double */
	double scl = 1.0*M_PI;
	cuchebCheckError(dfuncaller(dn, dpts, 1, dfvs, 1, (void*)&scl),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* double */
	cuchebCheckError(cudaMemcpy(ddfvs, dfvs, dn*sizeD, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebcoeffs
	/* double */
	cuchebCheckError(cuchebDcoeffs(dn, ddfvs, 1, ddcfs, 1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* double */
	cuchebCheckError(cudaMemcpy(dcfs, ddcfs, dn*sizeD, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* double */
	printf("\ndouble precision\n");
	for(int ii=0;ii<dn;ii++){printf("dpts[%d] = %+1.15e, dfvs[%d] = %+1.15e, dcfs[%d] = %+1.15e\n",ii,dpts[ii],ii,dfvs[ii],ii,dcfs[ii]);}
	printf("\n");

	// free device memory
	/* double */
	cuchebCheckError(cudaFree(ddcfs),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ddfvs),__FILE__,__LINE__);

	// free host memory
	/* double */
	free(dpts);
	free(dcfs);
	free(dfvs);
	
	return 0;
}
