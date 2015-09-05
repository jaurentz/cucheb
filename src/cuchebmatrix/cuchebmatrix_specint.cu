#include <cucheb.h>

/* routine to estimate spectral interval */
int cuchebmatrix_specint(cuchebmatrix* ccm){

  // number of arnoldi steps
  int nvecs;
  nvecs = min(ccm->m,200);

  // create cuchebpoly object
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // create lanczos object
  cucheblanczos ccl;
  cucheblanczos_init(nvecs,ccm,&ccl);

  // set starting vector
  cucheblanczos_startvec(&ccl);

  // arnoldi run
  cucheblanczos_arnoldi(ccm,&ccl);

  // compute ritz values
  cucheblanczos_eig(ccm,&ccl);

  // compute residuals
  double a, b;
  double ar, br;
  a = (ccl.diag)[nvecs-1];
  b = (ccl.diag)[0];
  ar = (ccl.sdiag)[nvecs-1]*abs(ccl.schurvecs[nvecs*nvecs-1]);
  br = (ccl.sdiag)[nvecs-1]*abs(ccl.schurvecs[nvecs-1]);
  ccm->a = a;
  ccm->b = b;
//  printf(" a = %+e, ar = %+e\n", a, ar);
//  printf(" b = %+e, br = %+e\n", b, br);

  // refine lower bound if necessary
  double tol = 1e-4;
  if ( ar >= tol*max(abs(a),abs(b)) ) {
  
    // create point filter
    cuchebpoly_pointfilter(a,b,a,50,&ccp);

    // set starting vector
    cucheblanczos_startvec(&ccl);

    // arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,&ccl);

    // compute ritz values
    cucheblanczos_eig(ccm,&ccl);

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,&ccl);

    // compute new candidate b
    a = (ccl.diag)[0];
    ar = (ccl.sdiag)[0];
    ccm->a = a;
    printf(" a = %+e, ar = %+e\n", a, ar);

  }

  // refine upper bound if necessary
  if ( br >= tol*max(abs(a),abs(b)) ) {
  
    // create point filter
    cuchebpoly_pointfilter(a,b,b,50,&ccp);

    // set starting vector
    cucheblanczos_startvec(&ccl);

    // arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,&ccl);

    // compute ritz values
    cucheblanczos_eig(ccm,&ccl);

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,&ccl);

    // compute new candidate b
    b = (ccl.diag)[0];
    br = (ccl.sdiag)[0];
    ccm->b = b;
    printf(" b = %+e, br = %+e\n", b, br);

  }

  // add a litte wiggle room
//  double eps;
//  eps = tol*max(abs(ccm->a),abs(ccm->b));
//  ccm->a = ccm->a - eps;
//  ccm->b = ccm->b + eps;

  // destroy ccl
  cucheblanczos_destroy(&ccl);

  // destory CCP
  cuchebpoly_destroy(&ccp);

  // return  
  return 0;

}
