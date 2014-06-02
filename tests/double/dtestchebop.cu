#include <cucheb.h>
#include <omp.h>

/* helpers for double precision constructors */
__global__ void dfunkernel(int n, const double *in, int incin, double* out, int incout, void* userdata);
cuchebStatus_t dfuncaller(int n, const double *in, int incin, double* out, int incout, void* userdata);
__global__ void testopkernel(int n, double *x, double *y, double a, double b);
void testop(void *x, void *y, void *user);

class Lap{
	public:
		int n;
		double a;
		double b;
};

// driver
int main(void){
	
	// compute variables
	Lap LD;
	double *x, *y, *dx, *dy;
	
	// begin timer
	double begin, end;
	begin = omp_get_wtime();

	// set LD
	LD.n = pow(2,5);
	LD.a = 0.0;
	LD.b = 2.0*(double)(LD.n+1);
	printf("\nn = %d\n",LD.n);
	
	// set chebpoly
	double tol = 1e-2;
	void* userdata;
	ChebPoly CP(&dfuncaller,&LD.a,&LD.b,userdata,&tol);
	CP.print();

	// allocate memory
	cuchebCheckError((void*)(x = (double*)malloc((LD.n)*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError((void*)(y = (double*)malloc(LD.n*sizeof(double))),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dx,LD.n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dy,LD.n*sizeof(double)),__FILE__,__LINE__);
	
	// initialize memory
	cuchebCheckError(cuchebDinit(LD.n,dx,1,1.0),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(LD.n,dy,1,0.0),__FILE__,__LINE__);
	
	// multiply
	cuchebCheckError(cuchebDmult(LD.n,dx,dy,&testop,(void*)&LD,&CP),__FILE__,__LINE__);
	
	// copy memory to host
	cuchebCheckError(cudaMemcpy(x,dx,LD.n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(y,dy,LD.n*sizeof(double),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
	// print
	for(int ii=0;ii<10;ii++){
		printf("x[%d] = %+e, y[%d] = %+e\n",ii,x[ii],ii,y[ii]);
	}
	
	// end timer
	end = omp_get_wtime();
	printf("\ntime = %e\n\n",end-begin);

	// free memory
	free(x);
	free(y);
	cuchebCheckError(cudaFree(dx),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(dy),__FILE__,__LINE__);

	return 0;
}


/* helpers for single precision constructors */
/* kernel to call devince function */
__global__ void dfunkernel(int n, const double *in, int incin, double* out, int incout, void* userdata){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

	if(ii < n){
		out[ii*incout] = exp(-(1e2)*in[ii*incin]);
		//out[ii*incout] = 1.0f*in[ii*incin]*in[ii*incin];
	}
}
/* subroutine to call sfunkernel */
cuchebStatus_t dfuncaller(int n, const double *in, int incin, double* out, int incout, void* userdata){
	
	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check incin
	if(incin <= 0){
		fprintf(stderr,"\nIn %s line: %d, incin must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check incout
	if(incout <= 0){
		fprintf(stderr,"\nIn %s line: %d, incout must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// launch fill input kernel
	dfunkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout, userdata);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}

__global__ void testopkernel(int n, double *x, double *y, double a, double b){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	double scl = (b-a)/(n+1);

	if(ii < n){
		//if(ii == 0){y[ii] = (2.0*x[ii] - x[ii+1])/scl/scl;}
		//else if(ii == n-1){y[ii] = (2.0*x[ii] - x[ii-1])/scl/scl;}
		//else{y[ii] = (2.0*x[ii] - x[ii+1] - x[ii-1])/scl/scl;}
		y[ii] = (a + (b-a)*((double)ii/(n-1)))*x[ii];
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

