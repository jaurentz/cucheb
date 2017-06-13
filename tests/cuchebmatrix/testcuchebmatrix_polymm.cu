#include <cucheb.h>

/* driver */
int main(){

  // input file
  string mtxfile("../matrices/SiH4.mtx");

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
  cuchebpoly_pointfilter(ccm.a,ccm.b,0.0,100,&ccp);
  cuchebpoly_print(&ccp);

  // create some vectors on the GPU
  int bsize = 3;
  double alpha, beta;
  double* dx;
  double* dy;
  double* dv1;
  double* dv2;

  cudaMalloc(&dx,bsize*(ccm.n)*sizeof(double));
  cudaMalloc(&dy,bsize*(ccm.n)*sizeof(double));
  cudaMalloc(&dv1,bsize*(ccm.n)*sizeof(double));
  cudaMalloc(&dv2,bsize*(ccm.n)*sizeof(double));

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

  // compute dy = p(A)*dx
  cuchebmatrix_polymm(&ccm,&ccp,bsize,dx,dy,dv1,dv2);

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

  // destroy CCM
  cuchebmatrix_destroy(&ccm);

  // destroy CCP
  cuchebpoly_destroy(&ccp);

  // destroy vectors
  cudaFree(dx);
  cudaFree(dy);
  cudaFree(dv1);
  cudaFree(dv2);

  // return 
  return 0;

}
