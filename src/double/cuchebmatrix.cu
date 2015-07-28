#include <cuchebmatrix.h>

/* routine to initialize cuchebmatrix object */
int cuchebmatrix_init(char* mtxfile, cuchebmatrix* ccm){

  // compute variables
  int ret;
  FILE *fptr;

  // attempt to open file
  fptr = fopen(mtxfile, "r");
  if (fptr == NULL) { 
    printf("Could not open matrix file.\n");
    exit(1); 
  }

  // attempt to read file banner
  ret = mm_read_banner(fptr, &(ccm->matcode));
  if (ret != 0) {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }

  // get matrix dimensions
  ret = mm_read_mtx_crd_size(fptr, &(ccm->m), &(ccm->n), &(ccm->nnz));
  if (ret != 0) {
    printf("Could not read matrix dimensions.\n");
    exit(1);
  }

  // allocate memory
  ccm->rowinds = new int[ccm->nnz];
  if (ccm->rowinds == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccm->colinds = new int[ccm->nnz];
  if (ccm->colinds == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccm->vals = new double[ccm->nnz];
  if (ccm->vals == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // read in matrix entries
  for (int ii=0; ii<ccm->nnz; ii++) {
    ret = fscanf(fptr, "%d %d %lg\n", &(ccm->rowinds)[ii], &(ccm->colinds)[ii],
                 &(ccm->vals)[ii]);
    (ccm->rowinds)[ii]--;  
    (ccm->colinds)[ii]--;
  }

  // close file
  if (fptr != stdin) fclose(fptr);

  // return  
  return 0;

}

/* routine to free memory in cuchebmatrix object */
int cuchebmatrix_destroy(cuchebmatrix* ccm){

  // free rowinds
  delete[] ccm->rowinds;

  // free colinds
  delete[] ccm->colinds;

  // free vals
  delete[] ccm->vals;
 
  // return  
  return 0;

}

/* routine for standard print */
int cuchebmatrix_print(cuchebmatrix* ccm){

  // print banner
  printf("\ncuchebmatrix:\n");

  // print matcode
  printf(" matcode = %s\n",mm_typecode_to_str(ccm->matcode));
 
  // print m
  printf(" m = %d\n",ccm->m);
 
  // print n
  printf(" n = %d\n",ccm->n);
 
  // print nnz
  printf(" nnz = %d\n",ccm->nnz);
  printf("\n");
 
  // return 
  return 0;

}

/* routine for long print */
int cuchebmatrix_printlong(cuchebmatrix* ccm){

  // print banner
  printf("\ncuchebmatrix:\n");

  // print matcode
  printf(" matcode = %s\n",mm_typecode_to_str(ccm->matcode));
 
  // print m
  printf(" m = %d\n",ccm->m);
 
  // print n
  printf(" n = %d\n",ccm->n);
 
  // print nnz
  printf(" nnz = %d\n",ccm->nnz);

  // print rowinds, colinds and vals
  for (int ii=0; ii<ccm->nnz; ii++) {
    printf(" rowinds[%d] = %d, colinds[%d] = %d, vals[%d] = %+e\n",
           ii,(ccm->rowinds)[ii],ii,(ccm->colinds)[ii],ii,(ccm->vals)[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

