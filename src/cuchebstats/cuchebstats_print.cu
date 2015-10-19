#include <cucheb.h>

/* routine for standard print */
int cuchebstats_print(cuchebstats* ccs){

  // print banner
  printf("\ncuchebstats:\n");

  // print members
  printf(" mat_dim = %d\n",ccs->mat_dim);
  printf(" mat_nnz = %d\n",ccs->mat_nnz);
  printf("\n");
  printf(" block_size     = %d\n",ccs->block_size);
  printf(" num_blocks     = %d\n",ccs->num_blocks);
  printf(" num_iters      = %d\n",ccs->num_iters);
  printf(" num_innerprods = %d\n",ccs->num_innerprods);
  printf(" max_degree     = %d\n",ccs->max_degree);
  printf(" num_matvecs    = %d\n",ccs->num_matvecs);
  printf(" specint_time   = %e (sec)\n",ccs->specint_time);
  printf(" innerprod_time   = %e (sec)\n",ccs->innerprod_time);
  printf(" matvec_time   = %e (sec)\n",ccs->matvec_time);
  printf("\n");
  printf(" num_conv = %d\n",ccs->num_conv);
  printf(" max_res  = %e\n",ccs->max_res);
  printf("\n");

  // return 
  return 0;

}

