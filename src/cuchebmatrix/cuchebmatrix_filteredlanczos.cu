#include <cucheb.h>

/* filtered lanczos routine for point value */
int cuchebmatrix_filteredlanczos(int neig, double shift, int bsize, cuchebmatrix* ccm,
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

  // initialize lanczos object
  cucheblanczos_init(bsize,MAX_NUM_BLOCKS,ccm,ccl);

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // initialize number of converged eigenvalues
  int numconv = 0;

  // loop through various filters
  double tau;
  tau = 10.0*(ccm->m);
  for (int jj=0; jj<MAX_RESTARTS+1; jj++) {

    // create filter polynomial
    cuchebpoly_gaussianfilter(ccm->a,ccm->b,rho,pow(10.0,jj)*tau,&ccp);
    cuchebpoly_print(&ccp);

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,ccl);

    // compute ritz values
    cucheblanczos_eig(ccm,ccl);

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,ccl);

    // check convergence
    cucheblanczos_checkconvergence(&numconv,rho,ccm,ccl); 

    // print eigenvalues
    for(int ii=0; ii < numconv; ii++){
      printf(" evals[%d] = %+e, res[%d] = %+e\n",
             ii,ccl->evals[ccl->index[ii]],ii,ccl->res[ccl->index[ii]]);
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





/* filtered lanczos routine for interval */
int cuchebmatrix_filteredlanczos(double lbnd, double ubnd, int bsize, cuchebmatrix* ccm, 
                                 cucheblanczos* ccl){

  // compute spectral interval
  cuchebmatrix_specint(ccm);
  cuchebmatrix_print(ccm);

  // make sure lbnd is valid
  if (isnan(lbnd)) {
    printf("lbnd cannot be NaN!\n");
    exit(1);
  }

  // make sure ubnd is valid
  if (isnan(ubnd)) {
    printf("ubnd cannot be NaN!\n");
    exit(1);
  }

  // check c and d
  if ( lbnd >= ubnd ) {
    return 1;
  }

  // compute lower bound 
  double a, b;
  a = ccm->a;
  b = ccm->b;
  double lb;
  if (lbnd <= a) {lb = a;}
  else if (lbnd >= b) {return 1;}
  else {lb = lbnd;}

  // compute upper bound 
  double ub;
  if (ubnd >= b) {ub = b;}
  else if (ubnd <= a) {return 1;}
  else {ub = ubnd;}

  // initialize lanczos object
  cucheblanczos_init(bsize,MAX_NUM_BLOCKS,ccm,ccl);

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // initialize number of converged eigenvalues
  int numconv = 0;

  // loop through various filters
  for (int jj=0; jj<MAX_RESTARTS+1; jj++) {

    // create filter polynomial
    cuchebpoly_stepfilter(ccm->a,ccm->b,lb,ub,50*(jj+2),&ccp);
    cuchebpoly_print(&ccp);

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,ccl);

    // compute ritz values
    cucheblanczos_eig(ccm,ccl);

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,ccl);

    // check convergence
    cucheblanczos_checkconvergence(&numconv,lb,ub,ccm,ccl); 

    // print eigenvalues
    for(int ii=0; ii < numconv; ii++){
      printf(" evals[%d] = %+e, res[%d] = %+e\n",
             ii,ccl->evals[ccl->index[ii]],ii,ccl->res[ccl->index[ii]]);
    }
    printf("\n");

    // exit if converged
    if (numconv > 0) { break; }

  }


  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // return  
  return 0;

}
