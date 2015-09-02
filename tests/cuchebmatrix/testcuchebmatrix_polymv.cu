#include <cucheb.h>

/* driver */
int main(){

  // input file
  //string mtxfile("./matrices/Sandi_authors.mtx");
  //string mtxfile("./matrices/Trefethen_20.mtx");
  string mtxfile("./matrices/Stranke94.mtx");

  // cuhebmatrix
  cuchebmatrix ccm;
  cuchebmatrix_init(mtxfile, &ccm);

  // compute spectral interval 
  cuchebmatrix_specint(&ccm);

  // print 
  cuchebmatrix_print(&ccm);

  // create a point filter
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);
  cuchebpoly_pointfilter(ccm.a,ccm.b,0.0,&ccp);
  cuchebpoly_print(&ccp);

  // create some vectors on the GPU
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

  // compute dy = p(A)*dx
  cuchebmatrix_polymv(&ccm,&ccp,dx,dy);

  // print dx and dy 
  for(int ii=0; ii < ccm.n; ii++){
    cudaMemcpy(&alpha,&dx[ii],sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(&beta,&dy[ii],sizeof(double),cudaMemcpyDeviceToHost);
    printf(" dx[%d] = %+e, dy[%d] = %+e\n", ii, alpha, ii, beta);
  }

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCP
  cuchebpoly_destroy(&ccp);

  // destroy vectors
  cudaFree(dx);
  cudaFree(dy);

  // return 
  return 0;

}
