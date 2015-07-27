#include <gpusollib.h>

__global__ void ddscal(const int n, const double alpha, double* y) {

    int i = blockIdx.x*blockDim.x + threadIdx.x;

    if( i < n ) y[i] = alpha *y [i];

}

