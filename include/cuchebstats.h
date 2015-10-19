#include <cuchebdependencies.h>

/* header file for cuchebstats data type */
#ifndef __cuchebstats_h__ 
#define __cuchebstats_h__

/* cuchebstats data type */
typedef struct {

  // problem dimensions
  int mat_dim;
  int mat_nnz;
  
  // iteration info
  int block_size;
  int num_blocks;
  int num_iters;
  int num_innerprods;
  int max_degree;
  int num_matvecs;
  double specint_time;
  double innerprod_time;
  double matvec_time;
  
  // convergence info
  int num_conv;
  double max_res;
 
} cuchebstats;

#endif /* __cuchebstats_h__ */
