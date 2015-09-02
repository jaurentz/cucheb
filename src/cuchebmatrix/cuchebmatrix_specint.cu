#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // number of arnoldi steps
  int nvecs;
  nvecs = min(ccm->m,20);

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
  a = (ccl.diag)[ccl.nvecs-1];
  b = (ccl.diag)[0];
  eps = .1*abs(b-a);
  ccm->a = a-eps;
  ccm->b = b+eps;

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // return  
  return 0;

}
