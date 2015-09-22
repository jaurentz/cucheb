#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // create stats object
  cuchebstats ccstats;

  // create lanczos object
  cucheblanczos ccl;
  cucheblanczos_init(1,DEF_NUM_VECS,ccm,&ccl);

  // set starting vector
  cucheblanczos_startvecs(&ccl);

  // loop to adaptively compute spectral interval
  int inda, indb;
  double nrm;
  int nres;
  nres = (DEF_NUM_VECS)/(DEF_STEP_SIZE) + 1;
  for (int jj=0; jj<nres; jj++){

    // arnoldi run
    cucheblanczos_arnoldi(DEF_STEP_SIZE,ccm,&ccl,&ccstats);

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

    // fudge factor
    ccm->a = ccm->a - pow(ccl.res[inda],2)*max(abs(ccm->a),abs(ccm->b));
    ccm->b = ccm->b + pow(ccl.res[indb],2)*max(abs(ccm->a),abs(ccm->b));

  }

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}
