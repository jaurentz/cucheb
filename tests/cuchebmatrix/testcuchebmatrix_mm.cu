#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("Trefethen_20.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;

  // initialize CCM
  cuchebmatrix_init(mtxfile, &ccm);

  // print 
  cuchebmatrix_print(&ccm);

  // create come vectors on the GPU
  int bsize = 3;
  double alpha, beta;
  double* dx;
  double* dy;

  cudaMalloc(&dx,bsize*(ccm.n)*sizeof(double));
  cudaMalloc(&dy,bsize*(ccm.n)*sizeof(double));

  // initialize dx to all 1's and dy to all zeros
  alpha = 1.0;
  beta = 0.0;
  for(int ii=0; ii < bsize*ccm.n; ii++){
    cudaMemcpy(&dx[ii],&alpha,sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(&dy[ii],&beta,sizeof(double),cudaMemcpyHostToDevice);
  }

  // print dx and dy 
  for(int ii=0; ii < ccm.n; ii++){
    for(int jj=0; jj < bsize; jj++){
      cudaMemcpy(&alpha,&dx[ccm.n*jj+ii],sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(&beta,&dy[ccm.n*jj+ii],sizeof(double),cudaMemcpyDeviceToHost);
      printf(" dx[%d] = %+e, dy[%d] = %+e, ", ccm.n*jj+ii, alpha, ccm.n*jj+ii, beta);
    }
      printf("\n");
  }
  printf("\n");

  // compute dy = alpha*ccm*dx + beta*dy
  alpha = 1.0;
  beta = 0.0;
  cuchebmatrix_mm(&ccm,bsize,&alpha,dx,&beta,dy);

  // print dx and dy 
  for(int ii=0; ii < ccm.n; ii++){
    for(int jj=0; jj < bsize; jj++){
      cudaMemcpy(&alpha,&dx[ccm.n*jj+ii],sizeof(double),cudaMemcpyDeviceToHost);
      cudaMemcpy(&beta,&dy[ccm.n*jj+ii],sizeof(double),cudaMemcpyDeviceToHost);
      printf(" dx[%d] = %+e, dy[%d] = %+e, ", ccm.n*jj+ii, alpha, ccm.n*jj+ii, beta);
    }
      printf("\n");
  }
  printf("\n");

  // compute specint
  cuchebmatrix_specint(&ccm);

  // print 
  cuchebmatrix_print(&ccm);

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy vectors
  cudaFree(dx);
  cudaFree(dy);

  // return 
  return 0;

}
