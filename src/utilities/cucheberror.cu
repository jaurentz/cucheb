/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>

/* cucheb exit */
void cuchebExit(int ii){
	cudaDeviceReset();
	exit(ii);
}

/* cucheb error string */
static const char* cuchebGetErrorString(cuchebStatus_t err){

	switch(err){
	    case CUCHEB_STATUS_SUCCESS:
	       return "CUCHEB_STATUS_SUCCESS"; 
	    case CUCHEB_STATUS_CUDA_FAILED:
	       return "CUCHEB_STATUS_CUDA_FAILED";
	    case CUCHEB_STATUS_CURAND_FAILED:
	       return "CUCHEB_STATUS_CURAND_FAILED"; 
	    case CUCHEB_STATUS_CUFFT_FAILED:
	       return "CUCHEB_STATUS_CUFFT_FAILED";
	    case CUCHEB_STATUS_CUBLAS_FAILED:
	       return "CUCHEB_STATUS_CUBLAS_FAILED"; 
	    case CUCHEB_STATUS_CUSPARSE_FAILED:
	       return "CUCHEB_STATUS_CUSPARSE_FAILED";
	    default: 
	       return "CUCHEB_STATUS_UNKNOWN";
	}
}

/* cublas error string */
static const char* cublasGetErrorString(cublasStatus_t err){

    switch (err){
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        default: 
	    	return "CUBLAS_STATUS_UNKNOWN";
    }
}

/* curand error string */
static const char* curandGetErrorString(curandStatus_t err)
{
    switch (err){
        case CURAND_STATUS_SUCCESS:
            return "CURAND_STATUS_SUCCESS";
        case CURAND_STATUS_VERSION_MISMATCH:
            return "CURAND_STATUS_VERSION_MISMATCH";
        case CURAND_STATUS_NOT_INITIALIZED:
            return "CURAND_STATUS_NOT_INITIALIZED";
        case CURAND_STATUS_ALLOCATION_FAILED:
            return "CURAND_STATUS_ALLOCATION_FAILED";
        case CURAND_STATUS_TYPE_ERROR:
            return "CURAND_STATUS_TYPE_ERROR";
        case CURAND_STATUS_OUT_OF_RANGE:
            return "CURAND_STATUS_OUT_OF_RANGE";
        case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
            return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
        case CURAND_STATUS_LAUNCH_FAILURE:
            return "CURAND_STATUS_LAUNCH_FAILURE";
        case CURAND_STATUS_PREEXISTING_FAILURE:
            return "CURAND_STATUS_PREEXISTING_FAILURE";
        case CURAND_STATUS_INITIALIZATION_FAILED:
            return "CURAND_STATUS_INITIALIZATION_FAILED";
        case CURAND_STATUS_ARCH_MISMATCH:
            return "CURAND_STATUS_ARCH_MISMATCH";
        case CURAND_STATUS_INTERNAL_ERROR:
            return "CURAND_STATUS_INTERNAL_ERROR";
        default:
            return "CURAND_STATUS_UNKNOWN";
    }
}

/* cufft error string */
static const char* cufftGetErrorString(cufftResult err)
{
    switch (err){
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            return "CUFFT_INCOMPLETE_PARAMETER_LIST";
        case CUFFT_INVALID_DEVICE:
            return "CUFFT_INVALID_DEVICE";
        case CUFFT_PARSE_ERROR:
            return "CUFFT_PARSE_ERROR";
        case CUFFT_NO_WORKSPACE:
            return "CUFFT_NO_WORKSPACE";
        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";
        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
        default: 
            return "CUFFT_UNKNOWN";
	}
}

/* cusparse error string */
static const char* cusparseGetErrorString(cusparseStatus_t err)
{
    switch (err){
        case CUSPARSE_STATUS_SUCCESS:
            return "CUSPARSE_STATUS_SUCCESS";
        case CUSPARSE_STATUS_NOT_INITIALIZED:
            return "CUSPARSE_STATUS_NOT_INITIALIZED";
        case CUSPARSE_STATUS_ALLOC_FAILED:
            return "CUSPARSE_STATUS_ALLOC_FAILED";
        case CUSPARSE_STATUS_INVALID_VALUE:
            return "CUSPARSE_STATUS_INVALID_VALUE";
        case CUSPARSE_STATUS_ARCH_MISMATCH:
            return "CUSPARSE_STATUS_ARCH_MISMATCH";
        case CUSPARSE_STATUS_MAPPING_ERROR:
            return "CUSPARSE_STATUS_MAPPING_ERROR";
        case CUSPARSE_STATUS_EXECUTION_FAILED:
            return "CUSPARSE_STATUS_EXECUTION_FAILED";
        case CUSPARSE_STATUS_INTERNAL_ERROR:
            return "CUSPARSE_STATUS_INTERNAL_ERROR";
        case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        default: 
            return "CUSPARSE_STATUS_UNKNOWN";
	}
}

/* host pointer */
void cuchebCheckError(void* err,char* file,int line){

	/* print error */
	if(err == NULL){
		fprintf(stderr,"\nHost memory allocation failure occured in %s at line: %d\n\n",file,line);
		cuchebExit(-1);
	}
} 

/* cucheb */
void cuchebCheckError(cuchebStatus_t err,char* file,int line){

	/* print error */
	if(err != CUCHEB_STATUS_SUCCESS){
		fprintf(stderr,"\n%s occured in %s at line: %d\n\n",cuchebGetErrorString(err),file,line);
		cuchebExit(-1);
	}
}  

/* cuda */
void cuchebCheckError(cudaError_t err,char* file,int line){

	/* print error and return CUDA_FAILED */
	if(err != cudaSuccess){
		fprintf(stderr,"\nCUCHEB_STATUS_CUDA_FAILED occured in %s at line: %d\n    cuda error: %s\n\n",
			file,line,cudaGetErrorString(err));
		cuchebExit(-1);
	}
} 

/* cublas */
void cuchebCheckError(cublasStatus_t err,char* file,int line){

	/* print error and return CUBLAS_FAILED */
	if(err != CUBLAS_STATUS_SUCCESS){
		fprintf(stderr,"\nCUCHEB_STATUS_CUBLAS_FAILED occured in %s at line: %d\n    cublas error: %s\n\n",
			file,line,cublasGetErrorString(err));
		cuchebExit(-1);
	}
}       

/* curand */
void cuchebCheckError(curandStatus_t err,char* file,int line){

	/* print error and return CURAND_FAILED */
	if(err != CURAND_STATUS_SUCCESS){
		fprintf(stderr,"\nCUCHEB_STATUS_CURAND_FAILED occured in %s at line: %d\n    curand error: %s\n\n",
			file,line,curandGetErrorString(err));
		cuchebExit(-1);
	}
} 

/* cufft */
void cuchebCheckError(cufftResult err,char* file,int line){

	/* print error and return CUFFT_FAILED */
	if(err != CUFFT_SUCCESS){
		fprintf(stderr,"\nCUCHEB_STATUS_CUFFT_FAILED occured in %s at line: %d\n    cufft error: %s\n\n",
			file,line,cufftGetErrorString(err));
		cuchebExit(-1);
	}
} 

/* cusparse */
void cuchebCheckError(cusparseStatus_t err,char* file,int line){

	/* print error and return CUSPARSE_FAILED */
	if(err != CUSPARSE_STATUS_SUCCESS){
		fprintf(stderr,"\nCUCHEB_STATUS_CUSPARSE_FAILED occured in %s at line: %d\n    cusparse error: %s\n\n",
			file,line,cusparseGetErrorString(err));
		cuchebExit(-1);
	}
} 
