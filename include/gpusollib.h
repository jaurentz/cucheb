#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cublas.h>
#include <assert.h>

// CUDA Runtime
#include <cuda_runtime.h>

// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cusparse_v2.h>
#include <cublas_v2.h>

#define IDX2C(i,j,ld) (((i)*(ld))+(j))
#define WARP      32
#define HALFWARP  16
#define BLOCKDIM  1024
#define ALIGNMENT 1024
#define MAXTHREADS (30 * 1024 * 5) // 1D kernels
// number of half warp per block, BLOCKDIM / 16
#define NHW 32
#define PI 3.14159265358979323846 // added by VK
#define ZERO 0.0
#define TRUE  1
#define FALSE 0
#define EPSILON   1.0e-18
#define EPSMAC    1.0e-16
/*--- max # of diags in DIA */
#define MAXDIAG 60
/*--- max # of color in 
 *--- multi-color reordering*/
#define MAXCOL 100
#include "datatype.h"
#include "protos.h"

