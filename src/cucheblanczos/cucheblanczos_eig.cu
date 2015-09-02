#include <cucheblanczos.h>

/* compute ritz values and vectors */
int cucheblanczos_eig(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int n, nvecs;
  double* diag;
  double* sdiag;
  double* schurvecs;
  double* dschurvecs;
  double* dvecs;
  n = ccl->n;
  nvecs = ccl->nvecs;
  diag = ccl->diag;
  sdiag = ccl->sdiag;
  schurvecs = ccl->schurvecs;
  dschurvecs = ccl->dschurvecs;
  dvecs = ccl->dvecs;

  // fill diag
  for(int ii=0; ii<nvecs; ii++){
    diag[ii] = schurvecs[ii*nvecs + ii];
  }

  // call lapack
  LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'i', nvecs, &diag[0], &sdiag[0],
                 &schurvecs[0], nvecs);

  // sort eigenvalues in descending order
  // create a vector of evals and indices
  vector< pair< double , int > > evals;
  for(int ii=0; ii < nvecs; ii++){
    evals.push_back(make_pair( -diag[ii] , ii ));
  }

  // sort vector
  sort(evals.begin(),evals.end());

  // update diag
  for(int ii=0; ii < nvecs; ii++){
    diag[ii] = -(evals[ii].first);
  }

  // update dschurvecs
  for(int ii=0; ii < nvecs; ii++){
    cudaMemcpy(&dschurvecs[ii*nvecs],&schurvecs[(evals[ii].second)*nvecs],
               nvecs*sizeof(double),cudaMemcpyHostToDevice);
  }

  // update schurvecs
  cudaMemcpy(&schurvecs[0],&dschurvecs[0],
               nvecs*nvecs*sizeof(double),cudaMemcpyDeviceToHost);

  // update dvecs
  double one = 1.0, zero = 0.0;
  cublasDgemm(ccm->cublashandle, CUBLAS_OP_N, CUBLAS_OP_N, n, nvecs, nvecs, &one,
              &dvecs[0], n, &dschurvecs[0], nvecs, &zero, &dvecs[0], n);

  // return  
  return 0;

}
