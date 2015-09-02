#include <cucheblanczos.h>

/* arnoldi run using cuchebmatrix */
int cucheblanczos_arnoldi(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, nvecs;
  double scl, one = 1.0, zero = 0.0, mone = -1.0;
  double* sdiag;
  double* dtemp;
  double* dvecs;
  double* dschurvecs;
  n = ccl->n;
  nvecs = ccl->nvecs;
  sdiag = ccl->sdiag;
  dtemp = ccl->dtemp;
  dvecs = ccl->dvecs;
  dschurvecs = ccl->dschurvecs;

  // loop through nvecs
  for(int ii=0; ii < nvecs; ii++){

    // apply matrix
    cuchebmatrix_mv(ccm,&one,&dvecs[ii*n],&zero,&dvecs[(ii+1)*n]);

    // orthogonalize
    cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, (ii+1), &one, &dvecs[0], n, 
                &dvecs[(ii+1)*n], 1, &zero, &dschurvecs[ii*nvecs], 1);
    cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, (ii+1), &mone, &dvecs[0], n, 
                &dschurvecs[ii*nvecs], 1, &one, &dvecs[(ii+1)*n], 1);

    // reorthogonalize
    cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, (ii+1), &one, &dvecs[0], n, 
                &dvecs[(ii+1)*n], 1, &zero, &dtemp[0], 1);
    cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, (ii+1), &mone, &dvecs[0], n, 
                &dtemp[0], 1, &one, &dvecs[(ii+1)*n], 1);
    cublasDaxpy(ccm->cublashandle, (ii+1), &one, &dtemp[0], 1, 
                &dschurvecs[ii*nvecs], 1);

    // normalize
    cublasDnrm2(ccm->cublashandle, n, &dvecs[(ii+1)*n], 1, &sdiag[ii]);
    scl = 1.0/sdiag[ii];
    cublasDscal(ccm->cublashandle, n, &scl, &dvecs[(ii+1)*n], 1);

  }

  // copy data to host
  cudaMemcpy(&(ccl->schurvecs)[0],&dschurvecs[0],nvecs*nvecs*sizeof(double),cudaMemcpyDeviceToHost);

  // return  
  return 0;

}

