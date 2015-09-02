#include <cucheb.h>

/* filtered lanczos routine */
int cuchebmatrix_filteredlanczos(int neig, double shift, cuchebmatrix* ccm,
                                 cucheblanczos* ccl){

  // compute spectral interval
  cuchebmatrix_specint(ccm);
  cuchebmatrix_print(ccm);

  // create filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);
  cuchebpoly_pointfilter(ccm->a,ccm->b,shift,&ccp);
  cuchebpoly_print(&ccp);

  // number of arnoldi steps
  int nvecs;
  nvecs = min(ccm->m,100);

  // initialize lanczos object
  cucheblanczos_init(nvecs,ccm,ccl);

  // set starting vector
  cucheblanczos_startvec(ccl);

  // filtered arnoldi run
  cucheblanczos_filteredarnoldi(ccm,&ccp,ccl);

  // compute ritz values
  cucheblanczos_eig(ccm,ccl);

  // compute rayleigh quotients
  cucheblanczos_rayleigh(ccm,ccl);

  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // return  
  return 0;

}
