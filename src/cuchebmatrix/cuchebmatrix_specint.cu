#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // create stats object
  cuchebstats ccstats;

  // create lanczos object
  cucheblanczos ccl;
  cucheblanczos_init(1,MAX_NUM_BLOCKS,ccm,&ccl);

  // set starting vector
  cucheblanczos_startvecs(&ccl);

  // loop to adaptively compute spectral interval
  int inda, indb;
  double nrm;
  for (int jj=0; jj<MAX_RESTARTS; jj++){

    // arnoldi run
    cucheblanczos_arnoldi(MAX_STEP_SIZE,ccm,&ccl,&ccstats);

    // compute ritz values
    cucheblanczos_ritz(ccm,&ccl);

    // set upper endpoint
    indb = ccl.index[0];
    ccm->b = ccl.evals[indb];

    // set lower endpoint
    inda = ccl.index[ccl.bsize*ccl.stop-1];
    ccm->a = ccl.evals[inda];

    // check convergence
    nrm = sqrt(DOUBLE_TOL)*max(abs(ccm->a),abs(ccm->b));
    if ( max(ccl.res[inda],ccl.res[indb]) < nrm ) { break; }

  }

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}
