#include <cucheb.h>

/* compute ritz values and vectors */
int cuchebblocklanczos_eig(cuchebmatrix* ccm, cuchebblocklanczos* ccb){

  // local variables
  int n, bsize, nblocks, nvecs;
  double* bands;
  double* evals;
  double* schurvecs;
  double* dschurvecs;
  double* dvecs;
  n = ccb->n;
  bsize = ccb->bsize;
  nblocks = ccb->nblocks;
  nvecs = bsize*nblocks;
  bands = ccb->bands;
  evals = ccb->evals;
  schurvecs = ccb->schurvecs;
  dschurvecs = ccb->dschurvecs;
  dvecs = ccb->dvecs;

  // fill bands
  for(int jj=0; jj<nvecs; jj++) {
    for(int ii=0; ii<bsize; ii++){
      bands[jj*(bsize+1)+ii] = schurvecs[jj*(nvecs+bsize)+jj+ii];
    }
  }

  // initialize schurvectors
  for(int jj=0; jj<nvecs; jj++) {
    for(int ii=0; ii<nvecs+bsize; ii++){
      if (ii == jj){ schurvecs[jj*(nvecs+bsize)+ii] = 1.0; }
      else{ schurvecs[jj*(nvecs+bsize)+ii] = 0.0; }
    }
  }

  // call lapack
  LAPACKE_dsbev(LAPACK_COL_MAJOR, 'v', 'l', nvecs, bsize, bands, bsize+1,
                evals, schurvecs, nvecs+bsize);

  // sort eigenvalues in descending order
  // create a vector of evals and indices
  vector< pair< double , int > > temp;
  for(int ii=0; ii < nvecs; ii++){
    temp.push_back(make_pair( -evals[ii] , ii ));
  }

  // sort vector
  sort(temp.begin(),temp.end());

  // update diag
  for(int ii=0; ii < nvecs; ii++){
    evals[ii] = -(temp[ii].first);
  }

  // update dschurvecs
  for(int ii=0; ii < nvecs; ii++){
    cudaMemcpy(&dschurvecs[ii*(nvecs+bsize)],&schurvecs[(temp[ii].second)*(nvecs+bsize)],
               (nvecs+bsize)*sizeof(double),cudaMemcpyHostToDevice);
  }

  // update schurvecs
  cudaMemcpy(&schurvecs[0],&dschurvecs[0],
               nvecs*(nvecs+bsize)*sizeof(double),cudaMemcpyDeviceToHost);

  // update dvecs
  double one = 1.0, zero = 0.0;
  cublasDgemm(ccm->cublashandle, CUBLAS_OP_N, CUBLAS_OP_N, n, nvecs, nvecs, &one,
              &dvecs[0], n, &dschurvecs[0], nvecs+bsize, &zero, &dvecs[0], n);

  // return  
  return 0;

}
