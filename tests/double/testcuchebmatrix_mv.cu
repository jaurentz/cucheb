#include <cuchebmatrix.h>

/* driver */
int main(){

  // input file
  //string mtxfile("./matrices/Sandi_authors.mtx");
  //string mtxfile("./matrices/Trefethen_20.mtx");
  string mtxfile("./matrices/Stranke94.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;

  // initialize CCM
  cuchebmatrix_init(mtxfile, &ccm);

  // print 
  cuchebmatrix_print(&ccm);

  // create come vectors on the GPU
  double alpha, beta;
  double* dx;
  double* dy;

  cudaMalloc(&dx,(ccm.n)*sizeof(double));
  cudaMalloc(&dy,(ccm.n)*sizeof(double));

  // initialize dx to all 1's and dy to all zeros
  alpha = 1.0;
  beta = 0.0;
  for(int ii=0; ii < ccm.n; ii++){
    cudaMemcpy(&dx[ii],&alpha,sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(&dy[ii],&beta,sizeof(double),cudaMemcpyHostToDevice);
  }

  // print dx and dy 
  for(int ii=0; ii < ccm.n; ii++){
    cudaMemcpy(&alpha,&dx[ii],sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(&beta,&dy[ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dx[%d] = %+e, dy[%d] = %+e\n", ii, alpha, ii, beta);
  }
  printf("\n");

  // compute dy = alpha*ccm*dx + beta*dy
  alpha = 1.0;
  beta = 0.0;
  cuchebmatrix_mv(&ccm,&alpha,dx,&beta,dy);

  // print dx and dy 
  for(int ii=0; ii < ccm.n; ii++){
    cudaMemcpy(&alpha,&dx[ii],sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(&beta,&dy[ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dx[%d] = %+e, dy[%d] = %+e\n", ii, alpha, ii, beta);
  }

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
