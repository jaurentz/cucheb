#include "gpusollib.h"

/*-----------------------------------------------*/
int cuda_init() {
  int deviceCount, dev;
  cudaDeviceProp deviceProp;
/*-----------------------------------------------*/
  cudaGetDeviceCount(&deviceCount);
  printf("=========================================\n");
  if (deviceCount == 0) {
    printf("There is no device supporting CUDA\n");
    return 1;
  }
  dev = 0;
  CUDA_SAFE_CALL(cudaSetDevice(dev));
  CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));
  if (deviceProp.major == 9999 && 
      deviceProp.minor == 9999) {
    printf("There is no device supporting CUDA.\n");
    return 1;
  }
  printf("Running on Device %d: \"%s\"\n", dev, deviceProp.name);
  printf("  Major revision number:          %d\n",
         deviceProp.major);
  printf("  Minor revision number:          %d\n",
           deviceProp.minor);
  printf("  Total amount of global memory:  %.2f GB\n",
         deviceProp.totalGlobalMem/1e9);
  printf("=========================================\n");
  return 0;
}

/*----------------------------------------*/
int gpusol_init() {
  if (cuda_init()) {
    printf("lol1");
    return 1;
  }

  if (cublasInit() != CUBLAS_STATUS_SUCCESS) {
    printf("lol2");
    return 2;
  }
  return 0;
}

/*--------------------------------------------------*/
void cuda_check_err() {
  cudaError_t cudaerr = cudaGetLastError() ;
  if (cudaerr != cudaSuccess) 
    printf("error: %s\n",cudaGetErrorString(cudaerr));
}

/*----------------------------------------------*/
int gpusol_finalize() {
/*------ Shut down CUBLAS */
  if (cublasShutdown() != CUBLAS_STATUS_SUCCESS) {
    printf("CUBLAS shut down FAILED !!\n");
    return 1;
  }
  cuda_check_err();
  return 0;
}

