#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // number of arnoldi steps
  int nblocks;
  nblocks = min(ccm->m,MAX_NUM_BLOCKS);

  // create lanczos object
  cucheblanczos ccl;
  cucheblanczos_init(1,nblocks,ccm,&ccl);

  // set starting vector
  cucheblanczos_startvecs(&ccl);

  // arnoldi run
  cucheblanczos_arnoldi(ccm,&ccl);

  // compute ritz values
  cucheblanczos_eig(ccm,&ccl);

  // compute residuals
  ccm->a = (ccl.evals)[nblocks-1];
  ccm->b = (ccl.evals)[0];

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}
