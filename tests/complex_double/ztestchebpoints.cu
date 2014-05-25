#include <cucheb.h>

int main(){

	// compute variables
	/* complex double */
	int zn = pow(2,3)+1;
	int sizeZ = sizeof(cuDoubleComplex);
	cuDoubleComplex *zpts, *za, *zb, *dzpts, *dza, *dzb;
	
	// allocate host memory
	/* complex double */
	cuchebCheckError((void*)(zpts = (cuDoubleComplex*)malloc(zn*sizeZ)),__FILE__,__LINE__);
	cuchebCheckError((void*)(za = (cuDoubleComplex*)malloc(sizeZ)),__FILE__,__LINE__);
	cuchebCheckError((void*)(zb = (cuDoubleComplex*)malloc(sizeZ)),__FILE__,__LINE__);
	
	// allocate device memory
	/* complex double */
	cuchebCheckError(cudaMalloc(&dzpts, zn*sizeZ),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dza, sizeZ),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dzb, sizeZ),__FILE__,__LINE__);
	
	// set host pointers
	/* complex double */
	*za = make_cuDoubleComplex(-1.0,0.0);
	*zb = make_cuDoubleComplex(1.0,0.0);
	
	// copy host memory to device memory
	/* complex double */
	cuchebCheckError(cudaMemcpy(dza, za, sizeZ, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(dzb, zb, sizeZ, cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// compute chebpoints
	/* complex double */
	cuchebCheckError(cuchebZpoints(zn,dza,dzb,dzpts,1),__FILE__,__LINE__);
	
	// copy device memory to host memory
	/* complex double */
	cuchebCheckError(cudaMemcpy(zpts, dzpts, zn*sizeZ, cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print output
	/* complex double */
	printf("complex double precision\n");
	for(int ii=0;ii<zn;ii++){printf("zpts[%d] = (%+1.15e,%+1.15e)\n",ii,cuCreal(zpts[ii]),cuCimag(zpts[ii]));}
	printf("\n");

	// free device memory
	/* complex double */
	cuchebCheckError(cudaFree(dzpts),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dza),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dzb),__FILE__,__LINE__);

	// free host memory
	/* complex double */
	free(zpts);
	free(za);
	free(zb);
	
	return 0;
}
