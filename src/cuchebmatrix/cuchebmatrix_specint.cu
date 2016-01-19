#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // create lanczos object
  cucheblanczos ccl;

  // call specint
  cuchebmatrix_specint(ccm,&ccl);

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm, cucheblanczos* ccl){

  // create stats object
  cuchebstats ccstats;

  // create lanczos object
  cucheblanczos_init(1,DEF_NUM_VECS,ccm,ccl);

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // loop to adaptively compute spectral interval
  int inda, indb;
  double nrm;
  int nres;
  nres = (DEF_NUM_VECS)/(DEF_STEP_SIZE) + 1;
  for (int jj=0; jj<nres; jj++){

    // arnoldi run
    cucheblanczos_arnoldi(DEF_STEP_SIZE,ccm,ccl,&ccstats);

    // compute ritz values
    cucheblanczos_ritz(ccm,ccl);

    // sort ritz values
    cucheblanczos_sort(ccl);

    // set upper endpoint
    indb = ccl->index[0];
    ccm->b = ccl->evals[indb];

    // set lower endpoint
    inda = ccl->index[ccl->bsize*ccl->stop-1];
    ccm->a = ccl->evals[inda];

    // fudge factor
    ccm->a = ccm->a - 10.0*sqrt(DOUBLE_TOL)*max(abs(ccm->a),abs(ccm->b));
    ccm->b = ccm->b + 10.0*sqrt(DOUBLE_TOL)*max(abs(ccm->a),abs(ccm->b));

    // check convergence
    nrm = sqrt(DOUBLE_TOL)*max(abs(ccm->a),abs(ccm->b));
    if ( max(ccl->res[inda],ccl->res[indb]) < nrm ) { break; }

  }

  // return  
  return 0;

}
