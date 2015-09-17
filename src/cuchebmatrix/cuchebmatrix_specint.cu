#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // number of arnoldi steps
  int nblocks;
  nblocks = min(ccm->m,200);

  // create lanczos object
  cucheblanczos ccl;
  cucheblanczos_init(1,nblocks,ccm,&ccl);

  // set starting vector
  cucheblanczos_startvecs(&ccl);

  // arnoldi run
  cucheblanczos_arnoldi(nblocks,ccm,&ccl);

  // compute ritz values
  cucheblanczos_ritz(ccm,&ccl);

  // set upper endpoint
  int indb = 0;
  ccm->b = -1e300;
  for (int ii=0; ii < ccl.nblocks; ii++){
    if (ccl.evals[ii] > ccm->b) {
      ccm->b = ccl.evals[ii];
      indb = ii;
    }
  }

  // set lower endpoint
  int inda = 0;
  ccm->a = ccm->b;
  for (int ii=0; ii < ccl.nblocks; ii++){
    if (ccl.evals[ii] < ccm->a) {
      ccm->a = ccl.evals[ii];
      inda = ii;
    }
  }

  // add fudge factor
  //ccm->a = ccm->a - ccl.res[inda]*max(abs(ccm->b),abs(ccm->a));
  //ccm->b = ccm->b + ccl.res[indb]*max(abs(ccm->b),abs(ccm->a));
  ccm->a = ccm->a - .01*abs(ccm->a);
  ccm->b = ccm->b + .01*abs(ccm->b);

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}
