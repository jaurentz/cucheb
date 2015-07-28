#include <cuchebmatrix.h>

/* routine to initialize cuchebmatrix object */
int cuchebmatrix_init(char* mtxfile, cuchebmatrix* ccm){

  // compute variables
  int ret;
  int new_nnz;
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

  // check matcode
  if ((ccm->matcode)[3] != 'S') {
    printf("Only symmetric matrices are supported.\n");
    exit(1);
  }

  // get matrix dimensions
  ret = mm_read_mtx_crd_size(fptr, &(ccm->m), &(ccm->n), &(ccm->nnz));
  if (ret != 0) {
    printf("Could not read matrix dimensions.\n");
    exit(1);
  }

  // allocate memory
  // for faster matvecs all elements of symmetric matrices must be stored
  ccm->rowinds = new int[2*(ccm->nnz)];
  if (ccm->rowinds == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccm->colinds = new int[2*(ccm->nnz)];
  if (ccm->colinds == NULL) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  ccm->vals = new double[2*(ccm->nnz)];
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

  // copy off diagonal elements
  new_nnz = 0;
  for (int ii=0; ii<ccm->nnz; ii++) {
    if ((ccm->rowinds)[ii] != (ccm->colinds)[ii]) {
      (ccm->rowinds)[(ccm->nnz)+new_nnz] = (ccm->colinds)[ii];
      (ccm->colinds)[(ccm->nnz)+new_nnz] = (ccm->rowinds)[ii];
      (ccm->vals)[(ccm->nnz)+new_nnz] = (ccm->vals)[ii];
      new_nnz++;
    }
  }
  ccm->nnz += new_nnz;

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

/* routine for sorting entries using GPU */
int cuchebmatrix_sort(cuchebmatrix* ccm){

  // device memory
  int* d_rowinds;
  int* d_colinds;
  double* d_vals;
  void *pBuffer = NULL; 
  int *P = NULL; 
  size_t pBufferSizeInBytes = 0; 

  // allocate device memory
  cudaMalloc(&d_rowinds, (ccm->nnz)*sizeof(int));
  cudaMalloc(&d_colinds, (ccm->nnz)*sizeof(int));
  cudaMalloc(&d_vals, (ccm->nnz)*sizeof(double));

  // copy memory to device
  cudaMemcpy(&d_rowinds, &(ccm->rowinds), (ccm->nnz)*sizeof(int),
             cudaMemcpyHostToDevice); 
  cudaMemcpy(&d_colinds, &(ccm->colinds), (ccm->nnz)*sizeof(int),
             cudaMemcpyHostToDevice); 
  cudaMemcpy(&d_vals, &(ccm->vals), (ccm->nnz)*sizeof(double),
             cudaMemcpyHostToDevice); 

  // Initialize CUSPARSE
  cusparseHandle_t cusparse_hand;
  cusparseCreate(&cusparse_hand);

  // step 1: allocate buffer
  cusparseXcoosort_bufferSizeExt(cusparse_hand, ccm->m, ccm->n, ccm->nnz,
                                 d_rowinds, d_colinds, &pBufferSizeInBytes); 
  cudaMalloc(&pBuffer, sizeof(char)* pBufferSizeInBytes);

// step 2: setup permutation vector P to identity 
//cudaMalloc( &P, sizeof(int)*nnz); 
//cusparseCreateIdentityPermutation(handle, nnz, P); 

// step 3: sort COO format by Row 
//cusparseXcoosortByRow(handle, m, n, nnz, cooRows, cooCols, P, pBuffer); 

// step 4: gather sorted cooVals 
//cusparseDgthr(handle, nnz, cooVals, cooVals_sorted, P, CUSPARSE_INDEX_BASE_ZERO);

  // Shutdown CUSPARSE
  cusparseDestroy(cusparse_hand);

  // copy memory to host
  cudaMemcpy(&(ccm->rowinds), &d_rowinds, (ccm->nnz)*sizeof(int),
             cudaMemcpyDeviceToHost); 
  cudaMemcpy(&(ccm->colinds), &d_colinds, (ccm->nnz)*sizeof(int),
             cudaMemcpyDeviceToHost); 
  cudaMemcpy(&(ccm->vals), &d_vals, (ccm->nnz)*sizeof(double),
             cudaMemcpyDeviceToHost); 

  // free device memory
  cudaFree(d_rowinds);
  cudaFree(d_colinds);
  cudaFree(d_vals);
  cudaFree(pBuffer);

  // return 
  return 0;

}

