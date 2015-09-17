#include <cucheb.h>

/* compute ritz values */
int cucheblanczos_ritz(cuchebmatrix* ccm, cucheblanczos* ccl){

  // local variables
  int bsize, nvecs, stop;
  int* index;
  double* evals;
  double* res;
  double* bands;
  double* schurvecs;
  bsize = ccl->bsize;
  nvecs = (ccl->bsize)*(ccl->nblocks);
  stop = ccl->stop;
  index = ccl->index;
  evals = ccl->evals;
  res = ccl->res;
  bands = ccl->bands;
  schurvecs = ccl->schurvecs;

  // reset index
  for(int ii=0; ii < stop*bsize; ii++){
    index[ii] = ii;
  }

  // fill bands
  for(int ii=0; ii<bsize+1; ii++){
    for(int jj=0; jj<stop*bsize; jj++) {
      bands[jj*(bsize+1)+ii] = schurvecs[jj*(nvecs+bsize)+jj+ii];
    }
  }

  // initialize residuals
  for(int ii=0; ii < stop*bsize; ii++){
    res[ii] = bands[(bsize+1)*stop*bsize-1];
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

  // update residuals
  for(int ii=0; ii < stop*bsize; ii++){
    res[ii] = res[ii]*abs(schurvecs[ii*(nvecs+bsize)+stop*bsize-1]);
  }

//  for(int ii=0; ii < stop*bsize; ii++){
//    printf(" %+e, %e\n",evals[ii],res[ii]);
//  }
//  printf("\n");

  // return  
  return 0;

}
