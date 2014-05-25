#include <cucheb.h>

int main(){

	// compute variables
	/* single */
	int sn = pow(2,3)+1;
	int sizeS = sizeof(float);
	float *spts, *sa, *sb, *dspts, *dsa, *dsb;
	
	// allocate host memory
	/* single */
	cuchebCheckError((void*)(spts = (float*)malloc(sn*sizeS)),__FILE__,__LINE__);
	cuchebCheckError((void*)(sa = (float*)malloc(sizeS)),__FILE__,__LINE__);
	cuchebCheckError((void*)(sb = (float*)malloc(sizeS)),__FILE__,__LINE__);
	
	// allocate device memory
	/* single */
	cuchebCheckError(cudaMalloc(&dspts, sn*sizeS),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dsa, sizeS),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dsb, sizeS),__FILE__,__LINE__);
	
	// set host pointers
	/* single */
	*sa = -1.0f;
	*sb = 1.0f;
	
	// copy host memory to device memory
	/* single */
	cuchebCheckError(cudaMemcpy(dsa, sa, sizeS, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dsb, sb, sizeS, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* single */
	cuchebCheckError(cuchebSpoints(sn,dsa,dsb,dspts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* single */
	cuchebCheckError(cudaMemcpy(spts, dspts, sn*sizeS, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* single */
	printf("single precision\n");
	for(int ii=0;ii<sn;ii++){printf("spts[%d] = %+e\n",ii,spts[ii]);}
	printf("\n");

	// free device memory
	/* single */
	cuchebCheckError(cudaFree(dspts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dsa),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dsb),__FILE__,__LINE__);

	// free host memory
	/* single */
	free(spts);
	free(sa);
	free(sb);
	
	return 0;
}
