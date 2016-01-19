#include <cucheb.h>

/* lanczos routine for interval */
int cuchebmatrix_lanczos(double lbnd, double ubnd, 
                         int bsize, int numvecs, int stepsize, 
                         cuchebmatrix* ccm, cucheblanczos* ccl){

  // cuchebstats variable
  cuchebstats ccstats;

  // call Lanczos
  cuchebmatrix_lanczos(lbnd,ubnd,bsize,numvecs,stepsize,ccm,ccl,&ccstats);

  // return
  return 0;

}


/* lanczos routine for interval with statistics */
int cuchebmatrix_lanczos(double lbnd, double ubnd, 
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
  ccstats->max_degree = 1;
  ccstats->num_matvecs = 0;
  ccstats->specint_time = 0.0;
  ccstats->innerprod_time = 0.0;
  ccstats->matvec_time = 0.0;
  ccstats->total_time = 0.0;
  ccstats->num_conv = 0;
  ccstats->max_res = 0.0;

  // collect some matrix statistics
  ccstats->mat_dim = ccm->m;
  ccstats->mat_nnz = ccm->nnz;

  // timing variables
  clock_t tick;

  // compute spectral interval
  tick = clock();
  cuchebmatrix_specint(ccm);

  // record compute time
  ccstats->specint_time = (clock()-tick);

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
    printf("\ncuchebmatrix_lanczos:\n");
    printf(" lb must be less than ub!\n\n");
    exit(1);
  }

  // initialize lanczos object
  cucheblanczos_init(bsize,numvecs,ccm,ccl);

  // collect some lanczos statistics
  ccstats->block_size = ccl->bsize;

  // set starting vector
  cucheblanczos_startvecs(ccl);

  // start stop watch
  tick = clock();

  // loop through various Krylov subspaces
  int numint = 0;
  int nres;
  int step;
  step = min(max(stepsize,1),ccl->nblocks);
  if ( (ccl->nblocks) % step == 0 ){ nres = (ccl->nblocks)/(step); }
  else{ nres = (ccl->nblocks)/(step) + 1; }

  for (int jj=0; jj<nres; jj++) {

    // filtered arnoldi run
    cucheblanczos_arnoldi(step,ccm,ccl,ccstats);

    // update ccstats
    // num_iters
    ccstats->num_iters += 1;

    // compute ritz values of A
    cucheblanczos_ritz(ccm,ccl);

    // sort evals
    cucheblanczos_sort(lb,ub,ccl);
    numint = ccl->nconv;

    // checkconvergence
    cucheblanczos_checkconvergence(ccl);

    // exit if converged
    if (ccl->nconv == numint && numint > 0) { break; }
    else { ccl->nconv = 0; }

  }

  // compute rayleigh quotients
  cucheblanczos_rayleigh(ccm,ccl);

  // num_conv
  ccstats->num_conv = ccl->nconv;

  // max_res
  for(int ii=0; ii < ccl->nconv; ii++){
    ccstats->max_res = max(ccstats->max_res,ccl->res[ccl->index[ii]]);
  }

  // record compute time
  ccstats->total_time = (clock()-tick) + ccstats->specint_time;
  ccstats->specint_time = (ccstats->specint_time)/((double)CLOCKS_PER_SEC);
  ccstats->total_time = (ccstats->total_time)/((double)CLOCKS_PER_SEC);

  // return  
  return 0;

}
