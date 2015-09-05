#include <cucheb.h>

/* filtered blocklanczos routine */
int cuchebmatrix_filteredblocklanczos(int neig, double shift, int bsize, cuchebmatrix* ccm,
                                 cuchebblocklanczos* ccb){

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

  // initialize blocklanczos object
  cuchebblocklanczos_init(bsize,MAX_NUM_BLOCKS,ccm,ccb);

  // set starting vector
  cuchebblocklanczos_startvecs(ccb);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // initialize number of converged eigenvalues
  int numconv = 0;

  // loop through various filters
  int nvecs;
  nvecs = (ccb->bsize)*(ccb->nblocks);
  double tau;
  tau = 10.0*(ccm->m);
  for (int jj=0; jj<MAX_RESTARTS+1; jj++) {
  //for (int jj=0; jj<1; jj++) {

    // create filter polynomial
    //cuchebpoly_pointfilter(ccm->a,ccm->b,rho,50*(jj+3),&ccp);
    cuchebpoly_gaussianfilter(ccm->a,ccm->b,rho,pow(10.0,jj)*tau,&ccp);
    cuchebpoly_print(&ccp);

    // filtered arnoldi run
    cuchebblocklanczos_filteredarnoldi(ccm,&ccp,ccb);

    // compute ritz values
    cuchebblocklanczos_eig(ccm,ccb);

    // print eigenvalues
//    for(int ii=0; ii < ccb->nvecs; ii++){
//      printf(" diag[%d] = %+e, sdiag[%d] = %+e\n",
//             ii,ccb->diag[ii],ii,ccb->sdiag[ii]);
//    }
//    printf("\n");

    // compute rayleigh quotients
    cuchebblocklanczos_rayleigh(ccm,ccb);

    // print eigenvalues
//    for(int ii=0; ii < nvecs; ii++){
//      printf(" evals[%d] = %+e, res[%d] = %+e\n",
//             ii,ccb->evals[ii],ii,ccb->res[ii]);
//    }
//    printf("\n");

    // check convergence
    cuchebblocklanczos_checkconvergence(&numconv,rho,ccm,ccb); 

    // print eigenvalues
    for(int ii=0; ii < numconv; ii++){
      printf(" evals[%d] = %+e, res[%d] = %+e\n",
             ii,ccb->evals[ccb->index[ii]],ii,ccb->res[ccb->index[ii]]);
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
