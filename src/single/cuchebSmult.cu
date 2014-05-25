#include <cucheb.h>

__global__ void sclenshaw(int n, float *coeff, float *a, float *b, float scl, float *vecin, float *temp2, float *vecout, float *temp1){
	int ii = (blockIdx.z*gridDim.y*gridDim.x + blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x*blockDim.y*blockDim.z 
			+ threadIdx.z*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
	float temp;
	float alpha, beta;

	if(ii < n){
		alpha = 2.0f/(*b-*a);
		beta = -(*b+*a)/(*b-*a);
		temp = vecout[ii];
		vecout[ii] = *coeff*vecin[ii] + scl*alpha*temp2[ii] + scl*beta*vecout[ii] - temp1[ii];
		temp1[ii] = temp;
	}
}

cuchebStatus_t cuchebSmult(int n,float* x,float* y,cuchebOpMult opmult,void* userdata,ChebPoly* chebpoly){

	// check n
	if(n <= 0){
		fprintf(stderr,"\nIn %s line: %d, n must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// check chebpoly field
	if(chebpoly->getField() != CUCHEB_FIELD_FLOAT){
		fprintf(stderr,"\nIn %s line: %d, field must be CUCHEB_FIELD_FLOAT.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// allocate work space
	float *temp1, *temp2;
	cuchebCheckError(cudaMalloc(&temp1, n*sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&temp2, n*sizeof(float)),__FILE__,__LINE__);
	
	// initialize workspace
	cuchebCheckError(cuchebSinit(n,temp1,1,0.0f),__FILE__,__LINE__);
	cuchebCheckError(cuchebSinit(n,temp2,1,0.0f),__FILE__,__LINE__);
	
	// get degree
	int m = chebpoly->getDegree();
	
	// get coeffs
	float *coeffs = (float*)chebpoly->getCoeffs();
	
	// get a and b
	float *a = (float*)chebpoly->getA();
	float *b = (float*)chebpoly->getB();
	
	// set blockSize and gridsize
	dim3 blockSize, gridSize;
	cuchebCheckError(cuchebSetGridBlocks(n,&blockSize,&gridSize),__FILE__,__LINE__);
	
	// set scl
	float scl = 1.0f;
	
	// degree 0 case
	if(m == 0){
		sclenshaw<<<gridSize,blockSize>>>(n, &coeffs[0], a, b, scl, x, temp2, y, temp1);
	}

	// arbitrary degree
	else if(m <= MAX_FLOAT_DEG){
		sclenshaw<<<gridSize,blockSize>>>(n, &coeffs[0], a, b, scl, x, temp2, y, temp1);
		cuchebCheckError(cuchebSinit(n,temp1,1,0.0f),__FILE__,__LINE__);

		scl = 2.0f;
		for(int jj = 0; jj < (m-1); jj++){
			opmult((void*)y,(void*)temp2,userdata);
			sclenshaw<<<gridSize,blockSize>>>(n, &coeffs[jj+1], a, b, scl, x, temp2, y, temp1);
		}

		scl = 1.0f;
		opmult((void*)y,(void*)temp2,userdata);
		sclenshaw<<<gridSize,blockSize>>>(n, &coeffs[m], a, b, scl, x, temp2, y, temp1);
	}
	
	// degree error
	else{
		fprintf(stderr,"\nIn %s line: %d, degree must be => 0 and <= %d.\n",__FILE__,__LINE__,MAX_FLOAT_DEG);
		cuchebExit(-1);
	}
	
	// free work space
	cuchebCheckError(cudaFree(temp1),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(temp2),__FILE__,__LINE__);
	
	// return
	return CUCHEB_STATUS_SUCCESS;
}
