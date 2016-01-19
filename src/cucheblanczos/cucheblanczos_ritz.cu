#include <cucheb.h>

/* compute ritz values */
int cucheblanczos_ritz(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int bsize, nvecs, stop;
  double V[MAX_BLOCK_SIZE] = {0.0};
  double R[(MAX_BLOCK_SIZE)*(MAX_BLOCK_SIZE)] = {0.0};
  double* evals;
  double* res;
  double* bands;
  double* schurvecs;
  bsize = ccl->bsize;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  stop = ccl->stop;
  evals = ccl->evals;
  res = ccl->res;
  bands = ccl->bands;
  schurvecs = ccl->schurvecs;

  // fill bands
  for(int ii=0; ii<bsize+1; ii++){
    for(int jj=0; jj<stop*bsize; jj++) {
      bands[jj*(bsize+1)+ii] = schurvecs[jj*(nvecs+bsize)+jj+ii];
    }
  }

  // fill R
  int ind;
  for(int ii=0; ii < bsize; ii++){
    ind = ((stop-1)*bsize+ii)*(nvecs+bsize+1)+bsize;
    for(int jj=0; jj < ii+1; jj++){
      R[ii*(MAX_BLOCK_SIZE + 1)-jj] = schurvecs[ind-jj];
    }
  }
  
  // initialize schurvectors
  for(int ii=0; ii<nvecs+bsize; ii++){
    for(int jj=0; jj<stop*bsize; jj++) {
      if (ii == jj){ schurvecs[jj*(nvecs+bsize)+ii] = 1.0; }
      else{ schurvecs[jj*(nvecs+bsize)+ii] = 0.0; }
    }
  }

  // call bandsymqr
  cuchebutils_bandsymqr(stop*bsize, bsize+1, bands, bsize+1,
                evals, schurvecs, nvecs+bsize);

  // compute residuals
  double tmp = 0;
  for(int ii=0; ii < stop*bsize; ii++){

    // set res[ii] to 0
    res[ii] = 0.0;

    // set V
    for(int jj=0; jj<bsize; jj++){
      V[jj] = schurvecs[ii*(nvecs+bsize) + (stop-1)*bsize + jj];
    }
  
    // compute matrix vector product and norm 
    for(int jj=0; jj<bsize; jj++){
      tmp = 0.0;
      for(int kk=0; kk<bsize-jj; kk++){
        tmp += R[kk*MAX_BLOCK_SIZE+jj]*V[kk];
      }
      res[ii] += tmp*tmp;
    }

    // compute sqrt
    res[ii] = sqrt(res[ii]);

  }

  // set nconv
  ccl->nconv = stop*bsize;

  // return  
  return 0;

}
