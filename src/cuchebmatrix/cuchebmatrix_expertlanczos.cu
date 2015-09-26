#include <cucheb.h>


/* expert lanczos routine for point with statistics */
int cuchebmatrix_expertlanczos(int neig, double shift, int degree,
                                 int bsize, int numvecs, int stepsize, 
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
    
  // initialize lanczos object
  cucheblanczos_init(bsize,numvecs,ccm,ccl);

  // collect some lanczos statistics
  ccstats->block_size = ccl->bsize;

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // compute interval
  double a, b, rho;
  a = ccm->a;
  b = ccm->b;
  rho = min(max(a,shift),b);
  double scl;
  scl = abs(b-a)*(10.0*neig*(ccl->bsize))/(ccm->m);
  double lb, ub;
  lb = max(a,rho-scl);
  ub = min(b,rho+scl);

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // create filter polynomial
  if (degree > -1) {
    cuchebpoly_stepfilter(ccm->a,ccm->b,lb,ub,degree,&ccp);
  }
  else {
    cuchebpoly_smartfilter(ccm->a,ccm->b,lb,ub,&ccp);
  }

  // max_degree
  ccstats->max_degree = max(ccstats->max_degree,ccp.degree);

  // start stop watch
  start = time(0);

  // loop through various Krylov subspaces
  int nres;
  int step;
  step = min(max(stepsize,1),ccl->nblocks);
  nres = (ccl->nblocks)/(step) + 1;
  for (int jj=0; jj<nres; jj++) {

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(step,ccm,&ccp,ccl,ccstats);

    // update ccstats
    // num_iters
    ccstats->num_iters += 1;

    // compute ritz values of p(A)
    cucheblanczos_ritz(ccm,ccl);

    // exit if converged
    if (ccl->nconv > neig) { 
      break; 
    }

  }

  // compute rayleigh quotients
  cucheblanczos_rayleigh(ccm,ccl);

  // sort evals
  cucheblanczos_sort(rho,ccl);

  // num_conv
  ccstats->num_conv = ccl->nconv;

  // max_res
  for(int ii=0; ii < ccl->nconv; ii++){
    ccstats->max_res = max(ccstats->max_res,ccl->res[ccl->index[ii]]);
  }

  // record compute time
  stop = time(0);
  ccstats->arnoldi_time = difftime(stop,start);

  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // return  
  return 0;

}




/* expert lanczos routine for interval with statistics */
int cuchebmatrix_expertlanczos(double lbnd, double ubnd, int degree,
                                 int bsize, int numvecs, int stepsize, 
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

  // make sure lbnd is valid
  if (isnan(lbnd)) {
    printf("lbnd cannot be NaN!\n");
    exit(1);
  }

  // compute lower bound 
  double a, b;
  a = ccm->a;
  b = ccm->b;
  double lb;
  lb = min(max(a,lbnd),b);

  // compute upper bound 
  double ub;
  ub = max(min(b,ubnd),a);

  // make sure ubnd is valid
  if (lb >= ub) {
    printf("\ncuchebmatrix_filteredlanczos:\n");
    printf(" lb must be less than ub!\n\n");
    exit(1);
  }

  // initialize filter polynomial
  cuchebpoly ccp;
  cuchebpoly_init(&ccp);

  // create filter polynomial
  if (degree > -1) {
    cuchebpoly_stepfilter(ccm->a,ccm->b,lb,ub,degree,&ccp);
  }
  else {
    cuchebpoly_smartfilter(ccm->a,ccm->b,lb,ub,&ccp);
  }

  // max_degree
  ccstats->max_degree = max(ccstats->max_degree,ccp.degree);
    
  // initialize lanczos object
  cucheblanczos_init(bsize,numvecs,ccm,ccl);

  // collect some lanczos statistics
  ccstats->block_size = ccl->bsize;

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // start stop watch
  start = time(0);

  // loop through various Krylov subspaces
  int numint = 0;
  int nres;
  int step;
  step = min(max(stepsize,1),ccl->nblocks);
  nres = (ccl->nblocks)/(step) + 1;
  for (int jj=0; jj<nres; jj++) {

    // filtered arnoldi run
    cucheblanczos_filteredarnoldi(step,ccm,&ccp,ccl,ccstats);

    // update ccstats
    // num_iters
    ccstats->num_iters += 1;

    // compute ritz values of p(A)
    cucheblanczos_ritz(ccm,ccl);

//  for(int ii=0; ii < ccl->nconv; ii++){
//    printf(" e,r = %+e,%e\n",ccl->evals[ccl->index[ii]],ccl->res[ccl->index[ii]]);
//  }
//printf("\n");

    // check to see if in interval
    numint = 0;
    for(int ii=0; ii<ccl->nconv; ii++){
      if(ccl->evals[ccl->index[ii]] >= .49){ numint += 1; }
      else { break; }
    }

    // exit if converged
    if (ccl->nconv > numint) { 
      ccl->nconv = numint;
      break; 
    }

  }

  // compute rayleigh quotients
  cucheblanczos_rayleigh(ccm,ccl);

  // sort evals
  cucheblanczos_sort(lb,ub,ccl);

  // num_conv
  ccstats->num_conv = ccl->nconv;

  // max_res
  for(int ii=0; ii < ccl->nconv; ii++){
    ccstats->max_res = max(ccstats->max_res,ccl->res[ccl->index[ii]]);
  }

  // record compute time
  stop = time(0);
  ccstats->arnoldi_time = difftime(stop,start);

  // destroy ccp
  cuchebpoly_destroy(&ccp);

  // return  
  return 0;

}
