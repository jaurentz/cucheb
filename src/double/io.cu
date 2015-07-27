#include <gpusollib.h>
//extern "C" {
#include <mmio.h>
//}
#define MAX_LINE 200

/*---------------------------------------------*
 *             READ COO Matrix Market          *
 *---------------------------------------------*/
int read_coo_MM(coo_t *coo, char *matfile) {

  FILE *p = fopen(matfile,"r");
  if (p == NULL) {
    printf("Unable to open file %s\n", matfile);
    exit(1);
  }
/* READ MM banner 
  if (mm_read_banner(p, &matcode) != 0){
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }
  if (!mm_is_valid(matcode)){
    printf("Invalid Matrix Market file.\n");
    exit(1);
  }
  if (!(mm_is_real(matcode) && mm_is_coordinate(matcode) 
        && mm_is_sparse(matcode))) {
    printf("Only sparse real-valued coordinate \
    matrices are supported\n");
    exit(1);
  }
*/
  int nrow, ncol, nnz, nnz2, k, j;
  char line[MAX_LINE];
/*------------- Read size */
  if (mm_read_mtx_crd_size(p, &nrow, &ncol, &nnz) !=0) {
    printf("MM read size error !\n");
    exit(1);
  }
  if (nrow != ncol) {
    fprintf(stdout,"This is not a square matrix!\n");
    exit(1);
  }

  nnz2 = nnz;
/*-------- Allocate mem for COO */
  coo->ir  = (int *)  malloc(nnz2 * sizeof(int));
  coo->jc  = (int *)  malloc(nnz2 * sizeof(int));
  coo->val = (double *) malloc(nnz2 * sizeof(double));
/*-------- read line by line */
  char *p1, *p2;
  for (k=0; k<nnz; k++) {
    fgets(line, MAX_LINE, p);
    for( p1 = line; ' ' == *p1; p1++ );
/*----------------- 1st entry - row index */
    for( p2 = p1; ' ' != *p2; p2++ ); 
    *p2 = '\0';
    float tmp1 = atof(p1);
    //coo->ir[k] = atoi(p1);
    coo->ir[k] = (int) tmp1;
/*-------------- 2nd entry - column index */
    for( p1 = p2+1; ' ' == *p1; p1++ );
    for( p2 = p1; ' ' != *p2; p2++ );
    *p2 = '\0';
    float tmp2 = atof(p1);
    coo->jc[k] = (int) tmp2;
    //coo->jc[k]  = atoi(p1);      
/*------------- 3rd entry - nonzero entry */
    p1 = p2+1;
    coo->val[k] = atof(p1); 
  }

  j = nnz;
  coo->n = nrow;
  coo->nnz = j;
  printf("matrix loaded ...\n");
  printf("  n = %d, nnz = %d\n", nrow, j);
  fclose(p);
  return 0;
}

