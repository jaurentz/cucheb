#include <cucheb.h>

int main(){

	// compute variables
	/* double */
	int dn = pow(2,3)+1;
	int sizeD = sizeof(double);
	double *dpts, *da, *db, *ddpts, *dda, *ddb;
	
	// allocate host memory
	/* double */
	cuchebCheckError((void*)(dpts = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(da = (double*)malloc(sizeD)),__FILE__,__LINE__);
	cuchebCheckError((void*)(db = (double*)malloc(sizeD)),__FILE__,__LINE__);
	
	// allocate device memory
	/* double */
	cuchebCheckError(cudaMalloc(&ddpts, dn*sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dda, sizeD),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&ddb, sizeD),__FILE__,__LINE__);
	
	// set host pointers
	/* double */
	*da = -1.0;
	*db = 1.0;
	
	// copy host memory to device memory
	/* double */
	cuchebCheckError(cudaMemcpy(dda, da, sizeD, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(ddb, db, sizeD, cudaMemcpyHostToDevice),__FILE__,__LINE__);

	// compute chebpoints
	/* double */
	cuchebCheckError(cuchebDpoints(dn,dda,ddb,ddpts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* double */
	cuchebCheckError(cudaMemcpy(dpts, ddpts, dn*sizeD, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* double */
	printf("double precision\n");
	for(int ii=0;ii<dn;ii++){printf("dpts[%d] = %+1.15e\n",ii,dpts[ii]);}
	printf("\n");

	// free device memory
	/* double */
	cuchebCheckError(cudaFree(ddpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dda),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(ddb),__FILE__,__LINE__);

	// free host memory
	/* double */
	free(dpts);
	free(da);
	free(db);
	
	return 0;
}
