#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // number of arnoldi steps
  int nvecs;
  nvecs = min(ccm->m,MAX_ARNOLDI_VECS);

  // create lanczos object
  cucheblanczos ccl;
  cucheblanczos_init(nvecs,ccm,&ccl);

  // set starting vector
  cucheblanczos_startvec(&ccl);

  // arnoldi run
  cucheblanczos_arnoldi(ccm,&ccl);

  // compute ritz values
  cucheblanczos_eig(ccm,&ccl);

  // estimate spectral interval
  double a, b, eps;
  double ar, br;
  a = (ccl.diag)[nvecs-1];
  b = (ccl.diag)[0];
  ar = (ccl.sdiag)[nvecs-1]*abs(ccl.schurvecs[nvecs*nvecs-1]);
  br = (ccl.sdiag)[nvecs-1]*abs(ccl.schurvecs[nvecs-1]);
  eps = .05*abs(b-a);
  ccm->a = a-eps;
  ccm->b = b+eps;

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}
