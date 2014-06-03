#include <cucheb.h>
#include <omp.h>

/* helpers for double precision constructors */
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
	double specrad;
	
	// begin timer
	double begin, end;
	begin = omp_get_wtime();

	// set LD
	LD.n = pow(2,20);
	LD.a = 0.0;
	LD.b = 1.0;//2.0*(double)(LD.n+1);
	printf("\nn = %d\n",LD.n);
	printf("\ntrue upperbound on specrad = %e\n",4.0*(LD.n+1.0)*(LD.n+1.0));

	// compute specrad
	cuchebCheckError(cuchebDspecrad(LD.n,&testop,(void*)&LD,&specrad),__FILE__,__LINE__);
	printf("\napproximate upperbound on specrad = %e\n",specrad);
	
	// end timer
	end = omp_get_wtime();
	printf("\ntime to compute specrad = %e\n\n",end-begin);

	return 0;
}


/* double precision ops */
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

