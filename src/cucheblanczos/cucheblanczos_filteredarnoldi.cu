#include <cucheb.h>

/* arnoldi run using cuchebmatrix */
int cucheblanczos_filteredarnoldi(int nsteps, cuchebmatrix* ccm, cuchebpoly* ccp,
                                  cucheblanczos* ccl, cuchebstats* ccstats){

  // local variables
  int n, bsize, nblocks, nvecs, stop;
  double scl, one = 1.0, zero = 0.0, mone = -1.0;
  double* dtemp;
  double* dv1;
  double* dv2;
  double* dschurvecs;
  double* dvecs;
  n = ccl->n;
  bsize = ccl->bsize;
  nblocks = ccl->nblocks;
  nvecs = bsize*nblocks;
  stop = ccl->stop;
  dtemp = ccl->dtemp;
  dschurvecs = ccl->dschurvecs;
  dvecs = ccl->dvecs;
  dv1 = ccl->dv1;
  dv2 = ccl->dv2;
  clock_t tick;

  // set niters
  int niters;
  niters = min(nsteps,nblocks-stop);

  // loop through nblocks
  int ind, odepth, start;
  for(int ii=0; ii < niters; ii++){

    // set index
    ind = (ii+stop)*bsize;

    // time matvecs
    tick = clock();

    // apply filtered matrix
    cuchebmatrix_polymm(ccm,ccp,bsize,&dvecs[ind*n],&dvecs[(ind+bsize)*n],dv1,dv2);
    cudaDeviceSynchronize();
    ccstats->matvec_time += (clock()-tick)/((double)CLOCKS_PER_SEC);

    // num_matvecs
    ccstats->num_matvecs += bsize*(ccp->degree);

    // inner loop for bsize blocks
    for(int jj=0; jj < bsize; jj++){

      // time innerprods
      tick = clock();

      // compute orthogonalization depth
      odepth = min((MAX_ORTH_DEPTH)*bsize+jj,ind+bsize);
      start = ind + jj + bsize - odepth;

      // orthogonalize
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, odepth, &one, &dvecs[start*n], 
                  n, &dvecs[(ind+bsize)*n], 1, &zero,
                  &dschurvecs[ind*(nvecs+bsize)+start], 1);
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, odepth, &mone, &dvecs[start*n], 
                  n, &dschurvecs[ind*(nvecs+bsize)+start], 1, &one, &dvecs[(ind+bsize)*n], 1);

      // num_innerprods 
      ccstats->num_innerprods += odepth;

      // reorthogonalize
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_T, n, odepth, &one, &dvecs[start*n], 
                  n, &dvecs[(ind+bsize)*n], 1, &zero, &dtemp[0], 1);
      cublasDgemv(ccm->cublashandle, CUBLAS_OP_N, n, odepth, &mone, &dvecs[start*n], 
                  n, &dtemp[0], 1, &one, &dvecs[(ind+bsize)*n], 1);
      cublasDaxpy(ccm->cublashandle, odepth, &one, &dtemp[0], 1, 
                  &dschurvecs[ind*(nvecs+bsize)+start], 1);

      // num_innerprods 
      ccstats->num_innerprods += odepth;

      // normalize
      cublasDnrm2(ccm->cublashandle, n, &dvecs[(ind+bsize)*n], 1, &scl);
      cudaMemcpy(&dschurvecs[ind*(nvecs+bsize)+ind+bsize], &scl,
                 sizeof(double), cudaMemcpyHostToDevice);
      scl = 1.0/scl;
      cublasDscal(ccm->cublashandle, n, &scl, &dvecs[(ind+bsize)*n], 1);
      cudaDeviceSynchronize();
      ccstats->innerprod_time += (clock()-tick)/((double)CLOCKS_PER_SEC);

    }

  }

  // update stop 
  ccl->stop += niters;

  // num_blocks
  ccstats->num_blocks += niters;

  // copy data to host
  cudaMemcpy(&(ccl->schurvecs)[0],&dschurvecs[0],
             (ccl->stop)*bsize*(nvecs+bsize)*sizeof(double),
             cudaMemcpyDeviceToHost);

  // return  
  return 0;

}

