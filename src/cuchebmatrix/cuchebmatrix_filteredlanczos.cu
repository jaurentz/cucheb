#include <cucheb.h>

/* filtered lanczos routine */
int cuchebmatrix_filteredlanczos(int neig, double shift, cuchebmatrix* ccm,
                                 cucheblanczos* ccl){

  // check neig
  if (neig > MAX_NUM_EIGS) {
    printf("Number of desired eigenvalues is too large!\n");
    exit(1);
  }

  // compute spectral interval
  cuchebmatrix_specint(ccm);
  cuchebmatrix_print(ccm);

  // make sure shift is valid
  double rho;
  if (isnan(shift)) {
    printf("Shift cannot be NaN!\n");
    exit(1);
  }
  else if (shift > ccm->b) { rho = ccm->b; }
  else if (shift < ccm->a) { rho = ccm->a; }
  else { rho = shift; }

  // number of arnoldi steps
  int nvecs;
  nvecs = min(ccm->m,MAX_ARNOLDI_VECS);

  // initialize lanczos object
  cucheblanczos_init(nvecs,ccm,ccl);

  // set starting vector
  cucheblanczos_startvec(ccl);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // initialize number of converged eigenvalues
  int numconv = 0;

  // loop through various filters
  double tau;
  tau = 10.0*(ccm->m);
  //for (int jj=0; jj<MAX_RESTARTS+1; jj++) {
  for (int jj=0; jj<1; jj++) {

    // create filter polynomial
    //cuchebpoly_pointfilter(ccm->a,ccm->b,shift,256,&ccp);
    cuchebpoly_gaussianfilter(ccm->a,ccm->b,rho,pow(10.0,jj)*tau,&ccp);
    cuchebpoly_print(&ccp);

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,ccl);

    // compute ritz values
    cucheblanczos_eig(ccm,ccl);

    // print eigenvalues
    for(int ii=0; ii < ccl->nvecs; ii++){
      printf(" diag[%d] = %+e, sdiag[%d] = %+e\n",
             ii,ccl->diag[ii],ii,ccl->sdiag[ii]);
    }
    printf("\n");

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,ccl);

    // print eigenvalues
    for(int ii=0; ii < ccl->nvecs; ii++){
      printf(" diag[%d] = %+e, sdiag[%d] = %+e\n",
             ii,ccl->diag[ii],ii,ccl->sdiag[ii]);
    }
    printf("\n");

    // check convergence
    cucheblanczos_checkconvergence(&numconv,rho,ccm,ccl); 

    // print eigenvalues
    for(int ii=0; ii < numconv; ii++){
      printf(" diag[%d] = %+e, sdiag[%d] = %+e\n",
             ii,ccl->diag[ccl->index[ii]],ii,ccl->sdiag[ccl->index[ii]]);
    }
    printf("\n");

    // exit if converged
    if (numconv >= neig) { break; }

  }


  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // return  
  return 0;

}
