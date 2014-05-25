#include <cucheb.h>

int main(){

	// compute variables
	/* complex single */
	int cn = pow(2,3)+1;
	int sizeC = sizeof(cuComplex);
	cuComplex *cpts, *ca, *cb, *dcpts, *dca, *dcb;
	
	// allocate host memory
	/* complex single */
	cuchebCheckError((void*)(cpts = (cuComplex*)malloc(cn*sizeC)),__FILE__,__LINE__);
	cuchebCheckError((void*)(ca = (cuComplex*)malloc(sizeC)),__FILE__,__LINE__);
	cuchebCheckError((void*)(cb = (cuComplex*)malloc(sizeC)),__FILE__,__LINE__);
	
	// allocate device memory
	/* complex single */
	cuchebCheckError(cudaMalloc(&dcpts, cn*sizeC),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dca, sizeC),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dcb, sizeC),__FILE__,__LINE__);
	
	// set host pointers
	/* complex single */
	*ca = make_cuFloatComplex(-1.0f,0.0f);
	*cb = make_cuFloatComplex(1.0f,0.0f);
	
	// copy host memory to device memory
	/* complex single */
	cuchebCheckError(cudaMemcpy(dca, ca, sizeC, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dcb, cb, sizeC, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* complex single */
	cuchebCheckError(cuchebCpoints(cn,dca,dcb,dcpts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex single */
	cuchebCheckError(cudaMemcpy(cpts, dcpts, cn*sizeC, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* complex single */
	printf("complex single precision\n");
	for(int ii=0;ii<cn;ii++){printf("cpts[%d] = (%+e,%+e)\n",ii,cuCrealf(cpts[ii]),cuCimagf(cpts[ii]));}
	printf("\n");

	// free device memory
	/* complex single */
	cuchebCheckError(cudaFree(dcpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dca),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dcb),__FILE__,__LINE__);
	
	// free host memory
	/* complex single */
	free(cpts);
	free(ca);
	free(cb);
	
	return 0;
}
