#include <cucheb.h>
#include <omp.h>

/* helpers for double precision constructors */
__global__ void testopkernel(int n, double *x, double *y, double a, double b);
void testop(void *x, void *y, void *user);

struct Lap{
		int n;
		double a;
		double b;
};

// driver
int main(void){
	
	// compute variables
	Lap LD;
	double *x, *b, *dx, *db;
	
	// begin timer
	double begin, end;
	begin = omp_get_wtime();

	// set LD
	LD.n = pow(2,5);
	LD.a = 0.0;
	LD.b = 2.0*(double)(LD.n+1);
	printf("\nn = %d\n",LD.n);

	// allocate memory
	cuchebCheckError((void*)(x = (double*)malloc((LD.n)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError((void*)(b = (double*)malloc(LD.n*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dx,LD.n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&db,LD.n*sizeof(double)),__FILE__,__LINE__);
	
	// initialize memory
	cuchebCheckError(cuchebDinit(LD.n,dx,1,0.0),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(LD.n,db,1,1.0),__FILE__,__LINE__);

	// solve
	cuchebDsolve(LD.n,&testop,(void*)&LD,dx,db);
	
	// copy memory to host
	cuchebCheckError(cudaMemcpy(x,dx,LD.n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(b,db,LD.n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print
	for(int ii=0;ii<LD.n;ii++){
		printf("x[%d] = %+e, b[%d] = %+e\n",ii,x[ii],ii,b[ii]);
	}
	
	// end timer
	end = omp_get_wtime();
	printf("\ntime = %e\n\n",end-begin);

	// free memory
	free(x);
	free(b);
	cuchebCheckError(cudaFree(dx),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(db),__FILE__,__LINE__);

	return 0;
}


/* helpers for double precision constructors */
__global__ void testopkernel(int n, double *x, double *y, double a, double b){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	double scl = (b-a)/(n+1);

	if(ii < n){
		if(ii == 0){y[ii] = (2.0*x[ii] - x[ii+1])/scl/scl;}
		else if(ii == n-1){y[ii] = (2.0*x[ii] - x[ii-1])/scl/scl;}
		else{y[ii] = (2.0*x[ii] - x[ii+1] - x[ii-1])/scl/scl;}
		//y[ii] = (a + (b-a)*((double)ii/(n-1)))*x[ii];
	}
}
void testop(void *x, void *y, void *user){
	
	Lap *LD = (Lap*)user;
	int n = LD->n;
	double a = LD->a;
	double b = LD->b;
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// launch fill input kernel
	testopkernel<<<gridSize,blockSize>>>(n, (double*)x, (double*)y, a, b);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
}

