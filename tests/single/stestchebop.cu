#include <cucheb.h>
#include <omp.h>

/* helpers for single precision constructors */
__global__ void sfunkernel(int n, const float *in, int incin, float* out, int incout);
cuchebStatus_t sfuncaller(int n, const float *in, int incin, float* out, int incout, void* userdata);
__global__ void testopkernel(int n, float *x, float *y, float a, float b);
void testop(void *x, void *y, void *user);

class Lap{
	public:
		int n;
		float a;
		float b;
};

// driver
int main(void){
	
	// compute variables
	Lap LD;
	float *x, *y, *dx, *dy;
	
	// begin timer
	double begin, end;
	begin = omp_get_wtime();

	// set LD
	LD.n = pow(2,20);
	LD.a = 0.0f;
	LD.b = 2.0f*(float)(LD.n+1);
	printf("\nn = %d\n",LD.n);
	
	// set chebpoly
	float tol = 1e-3;
	void * userdata;
	ChebPoly CP(&sfuncaller,&LD.a,&LD.b,userdata,&tol);
	CP.print();

	// allocate memory
	cuchebCheckError((void*)(x = (float*)malloc((LD.n)*sizeof(float))),__FILE__,__LINE__);
	cuchebCheckError((void*)(y = (float*)malloc(LD.n*sizeof(float))),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dx,LD.n*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&dy,LD.n*sizeof(float)),__FILE__,__LINE__);
	
	// initialize memory
	cuchebCheckError(cuchebSinit(LD.n,dx,1,1.0f),__FILE__,__LINE__);
	cuchebCheckError(cuchebSinit(LD.n,dy,1,0.0f),__FILE__,__LINE__);
	
	// multiply
	cuchebCheckError(cuchebSmult(LD.n,dx,dy,&testop,(void*)&LD,&CP),__FILE__,__LINE__);
	
	// copy memory to host
	cuchebCheckError(cudaMemcpy(x,dx,LD.n*sizeof(float),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	cuchebCheckError(cudaMemcpy(y,dy,LD.n*sizeof(float),cudaMemcpyDeviceToHost),__FILE__,__LINE__);
	
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
__global__ void sfunkernel(int n, const float *in, int incin, float* out, int incout){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;

	if(ii < n){
		out[ii*incout] = expf(-(1e3)*in[ii*incin]);
		//out[ii*incout] = 1.0f*in[ii*incin]*in[ii*incin];
	}
}
/* subroutine to call sfunkernel */
cuchebStatus_t sfuncaller(int n, const float *in, int incin, float* out, int incout, void* userdata){
	
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
	sfunkernel<<<gridSize,blockSize>>>(n, in, incin, out, incout);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}

__global__ void testopkernel(int n, float *x, float *y, float a, float b){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	float scl = (b-a)/(n+1);

	if(ii < n){
		if(ii == 0){y[ii] = (2.0f*x[ii] - x[ii+1])/scl/scl;}
		else if(ii == n-1){y[ii] = (2.0f*x[ii] - x[ii-1])/scl/scl;}
		else{y[ii] = (2.0f*x[ii] - x[ii+1] - x[ii-1])/scl/scl;}
		//y[ii] = (a + (b-a)*((float)ii/(n-1)))*x[ii];
	}
}
void testop(void *x, void *y, void *user){
	
	Lap *LD = (Lap*)user;
	int n = LD->n;
	float a = LD->a;
	float b = LD->b;
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// launch fill input kernel
	testopkernel<<<gridSize,blockSize>>>(n, (float*)x, (float*)y, a, b);
	
	// check for kernel error
	cuchebCheckError(cudaPeekAtLastError(),__FILE__,__LINE__);
}

