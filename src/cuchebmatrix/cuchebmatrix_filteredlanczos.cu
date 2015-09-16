#include <cucheb.h>

/* filtered lanczos routine for point value */
int cuchebmatrix_filteredlanczos(int neig, double shift, int bsize, cuchebmatrix* ccm,
                                 cucheblanczos* ccl){

  // temp variables
  cuchebstats ccstats;

  // call filtered lanczos
  return cuchebmatrix_filteredlanczos(neig,shift,bsize,ccm,ccl,&ccstats);

} 
  

/* filtered lanczos routine for point value */
int cuchebmatrix_filteredlanczos(int neig, double shift, int bsize, cuchebmatrix* ccm,
                                 cucheblanczos* ccl, cuchebstats* ccstats){

  // check neig
  if (neig > MAX_NUM_EIGS) {
    printf("\ncuchebmatrix_filteredlanczos:\n");
    printf(" Number of desired eigenvalues is too large!\n\n");
    exit(1);
  }

  // initialize ccstats
  ccstats->mat_dim = 0;
  ccstats->mat_nnz = 0;
  ccstats->block_size = 0;
  ccstats->num_blocks = 0;
  ccstats->num_iters = 0;
  ccstats->num_innerprods = 0;
  ccstats->max_degree = 0;
  ccstats->num_matvecs = 0;
  ccstats->specint_time = 0.0;
  ccstats->arnoldi_time = 0.0;
  ccstats->num_conv = 0;
  ccstats->max_res = 0.0;

  // collect some matrix statistics
  ccstats->mat_dim = ccm->m;
  ccstats->mat_nnz = ccm->nnz;

  // timing variables
  time_t start, stop;

  // compute spectral interval
  start = time(0);
  cuchebmatrix_specint(ccm);

  // record compute time
  stop = time(0);
  ccstats->specint_time = difftime(stop,start);

  // make sure shift is valid
  double rho;
  if (isnan(shift)) {
    printf("\ncuchebmatrix_filteredlanczos:\n");
    printf(" Shift cannot be NaN!\n\n");
    exit(1);
  }
  else if (shift > ccm->b) { rho = ccm->b; }
  else if (shift < ccm->a) { rho = ccm->a; }
  else { rho = shift; }

  // compute interval
  double lb, ub, scl;
  scl = .001*(ccm->b - ccm->a);
  lb = max(ccm->a,rho-scl);
  ub = min(ccm->b,rho+scl);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);
    
  // create filter polynomial
  cuchebpoly_stepfilter(ccm->a,ccm->b,lb,ub,50,&ccp);

  // initialize number of converged eigenvalues
  int numconv = 0;

  // start stop watch
  start = time(0);

  // loop through various filters
  for (int jj=0; jj<MAX_RESTARTS+1; jj++) {

    // initialize lanczos object
    cucheblanczos_init(bsize,(jj+2)*50,ccm,ccl);

    // collect some lanczos statistics
    ccstats->block_size = ccl->bsize;
    ccstats->num_blocks = ccl->nblocks;

    // set starting vector
    cucheblanczos_startvecs(ccl);

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,ccl,ccstats);

    // compute ritz values
    cucheblanczos_eig(ccm,ccl);

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,ccl);

    // check convergence
    cucheblanczos_checkconvergence(&numconv,rho,ccm,ccl); 

    // update ccstats
    // num_iters
    ccstats->num_iters += 1;

    // max_degree
    ccstats->max_degree = max(ccstats->max_degree,ccp.degree);
    
    // exit if converged
    if (numconv >= neig) { break; }

    // destroy ccl 
    if (jj < MAX_RESTARTS) {
      cucheblanczos_destroy(ccl);
    }

  }

  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // num_conv
  ccstats->num_conv = numconv;

  // max_res
  for(int ii=0; ii < numconv; ii++){
    ccstats->max_res = max(ccstats->max_res,ccl->res[ccl->index[ii]]);
  }
  ccstats->max_res = (ccstats->max_res)/max(abs(ccm->a),abs(ccm->b));

  // record compute time
  stop = time(0);
  ccstats->arnoldi_time = difftime(stop,start);

  // return  
  return 0;

}




/* filtered lanczos routine for interval */
int cuchebmatrix_filteredlanczos(double lbnd, double ubnd, int bsize, 
                                 cuchebmatrix* ccm, cucheblanczos* ccl){

  // temp variables
  cuchebstats ccstats;

  // call filtered lanczos
  return cuchebmatrix_filteredlanczos(lbnd,ubnd,bsize,ccm,ccl,&ccstats);

} 

/* filtered lanczos routine for interval with statistics */
int cuchebmatrix_filteredlanczos(double lbnd, double ubnd, int bsize, 
                                 cuchebmatrix* ccm, cucheblanczos* ccl, 
                                 cuchebstats* ccstats){

  // initialize ccstats
  ccstats->mat_dim = 0;
  ccstats->mat_nnz = 0;
  ccstats->block_size = 0;
  ccstats->num_blocks = 0;
  ccstats->num_iters = 0;
  ccstats->num_innerprods = 0;
  ccstats->max_degree = 0;
  ccstats->num_matvecs = 0;
  ccstats->specint_time = 0.0;
  ccstats->arnoldi_time = 0.0;
  ccstats->num_conv = 0;
  ccstats->max_res = 0.0;

  // collect some matrix statistics
  ccstats->mat_dim = ccm->m;
  ccstats->mat_nnz = ccm->nnz;

  // timing variables
  time_t start, stop;

  // compute spectral interval
  start = time(0);
  cuchebmatrix_specint(ccm);

  // record compute time
  stop = time(0);
  ccstats->specint_time = difftime(stop,start);
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

  // collect some lanczos statistics
  ccstats->block_size = ccl->bsize;
  ccstats->num_blocks = ccl->nblocks;

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // initialize number of converged eigenvalues
  int numconv = 0;

  // start stop watch
  start = time(0);

  // loop through various filters
  for (int jj=0; jj<MAX_RESTARTS+1; jj++) {
  //for (int jj=0; jj<1; jj++) {

    // create filter polynomial
    cuchebpoly_stepfilter(ccm->a,ccm->b,lb,ub,50*(jj+1),&ccp);

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(ccm,&ccp,ccl,ccstats);

    // compute ritz values
    cucheblanczos_eig(ccm,ccl);

    // compute rayleigh quotients
    cucheblanczos_rayleigh(ccm,ccl);

    // check convergence
    cucheblanczos_checkconvergence(&numconv,lb,ub,ccm,ccl); 

    // update ccstats
    // num_iters
    ccstats->num_iters += 1;

    // max_degree
    ccstats->max_degree = max(ccstats->max_degree,ccp.degree);
    
    // exit if converged
    if (numconv > 0) { break; }

  }

  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // num_conv
  ccstats->num_conv = numconv;

  // max_res
  for(int ii=0; ii < numconv; ii++){
    ccstats->max_res = max(ccstats->max_res,ccl->res[ccl->index[ii]]);
  }
  ccstats->max_res = (ccstats->max_res)/max(abs(ccm->a),abs(ccm->b));

  // record compute time
  stop = time(0);
  ccstats->arnoldi_time = difftime(stop,start);

  // return  
  return 0;

}
