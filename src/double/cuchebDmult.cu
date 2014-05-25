#include <cucheb.h>

__global__ void dclenshaw(int n, double *coeff, double *a, double *b, double scl, double *vecin, double *temp2, double *vecout, double *temp1){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	double temp;
	double alpha, beta;

	if(ii < n){
		alpha = 2.0/(*b-*a);
		beta = -(*b+*a)/(*b-*a);
		temp = vecout[ii];
		vecout[ii] = *coeff*vecin[ii] + scl*alpha*temp2[ii] + scl*beta*vecout[ii] - temp1[ii];
		temp1[ii] = temp;
	}
}

cuchebStatus_t cuchebDmult(int n,double* x,double* y,cuchebOpMult opmult,void* userdata,ChebPoly* chebpoly){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check chebpoly field
	if(chebpoly->getField() != CUCHEB_FIELD_DOUBLE){
		fprintf(stderr,"\nIn %s line: %d, field must be CUCHEB_FIELD_DOUBLE.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// allocate work space
	double *temp1, *temp2;
	cuchebCheckError(cudaMalloc(&temp1, n*sizeof(double)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&temp2, n*sizeof(double)),__FILE__,__LINE__);
	
	// initialize workspace
	cuchebCheckError(cuchebDinit(n,temp1,1,0.0),__FILE__,__LINE__);
	cuchebCheckError(cuchebDinit(n,temp2,1,0.0),__FILE__,__LINE__);
	
	// get degree
	int m = chebpoly->getDegree();
	
	// get coeffs
	double *coeffs = (double*)chebpoly->getCoeffs();
	
	// get a and b
	double *a = (double*)chebpoly->getA();
	double *b = (double*)chebpoly->getB();
	
	// initialize output
	cuchebCheckError(cuchebDinit(n,y,1,0.0),__FILE__,__LINE__);
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// set scl
	double scl = 1.0;
	
	// degree 0 case
	if(m == 0){
		dclenshaw<<<gridSize,blockSize>>>(n, &coeffs[0], a, b, scl, x, temp2, y, temp1);
	}

	// arbitrary degree
	else if(m <= MAX_DOUBLE_DEG){
		dclenshaw<<<gridSize,blockSize>>>(n, &coeffs[0], a, b, scl, x, temp2, y, temp1);
		cuchebCheckError(cuchebDinit(n,temp1,1,0.0),__FILE__,__LINE__);

		scl = 2.0;
		for(int jj = 0; jj < (m-1); jj++){
			opmult((void*)y,(void*)temp2,userdata);
			dclenshaw<<<gridSize,blockSize>>>(n, &coeffs[jj+1], a, b, scl, x, temp2, y, temp1);
		}

		scl = 1.0;
		opmult((void*)y,(void*)temp2,userdata);
		dclenshaw<<<gridSize,blockSize>>>(n, &coeffs[m], a, b, scl, x, temp2, y, temp1);
	}
	
	// degree error
	else{
		fprintf(stderr,"\nIn %s line: %d, degree must be => 0 and <= %d.\n",__FILE__,__LINE__,MAX_DOUBLE_DEG);
		cuchebExit(-1);
	}
	
	// free work space
	cuchebCheckError(cudaFree(temp1),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(temp2),__FILE__,__LINE__);
	
	// return
	return CUCHEB_STATUS_SUCCESS;
}
