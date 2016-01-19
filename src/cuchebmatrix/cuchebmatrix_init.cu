#include <cucheb.h>

/* routine to initialize cuchebmatrix object */
int cuchebmatrix_init(const string& mtxfile, cuchebmatrix* ccm){

  // compute variables
  int new_nnz;
  ifstream input_file;

  // attempt to open file
  input_file.open(mtxfile.c_str());
  if (!input_file.is_open()) { 
    printf("Could not open matrix file.\n");
    exit(1); 
  }

  // variables to parse banner
  string line, banner, mtx, crd, data_type, storage_scheme;

  // scan first line
  getline(input_file,line);
  istringstream iss(line);
  iss >> banner >> mtx >> crd >> data_type >> storage_scheme;

  // check for banner 
  if (banner.compare(MatrixMarketBanner) != 0)
    return MM_NO_HEADER;

  // first field should be "mtx"
  if (mtx.compare(MM_MTX_STR) != 0)
    return  MM_UNSUPPORTED_TYPE;
  mm_set_matrix(&(ccm->matcode));

  // second field describes whether this is a sparse matrix (in coordinate
  //  storage) or a dense array
  if (crd.compare(MM_SPARSE_STR) == 0)
    mm_set_sparse(&(ccm->matcode));
  else
  if (crd.compare(MM_DENSE_STR) == 0)
    mm_set_dense(&(ccm->matcode));
  else
    return MM_UNSUPPORTED_TYPE;
    
  // third field
  if (data_type.compare(MM_REAL_STR) == 0)
    mm_set_real(&(ccm->matcode));
  else
  if (data_type.compare(MM_COMPLEX_STR) == 0)
    mm_set_complex(&(ccm->matcode));
  else
  if (data_type.compare(MM_PATTERN_STR) == 0)
    mm_set_pattern(&(ccm->matcode));
  else
  if (data_type.compare(MM_INT_STR) == 0)
    mm_set_integer(&(ccm->matcode));
  else
    return MM_UNSUPPORTED_TYPE;
    
  // fourth field 
  if (storage_scheme.compare(MM_GENERAL_STR) == 0)
    mm_set_general(&(ccm->matcode));
  else
  if (storage_scheme.compare(MM_SYMM_STR) == 0)
    mm_set_symmetric(&(ccm->matcode));
  else
  if (storage_scheme.compare(MM_HERM_STR) == 0)
    mm_set_hermitian(&(ccm->matcode));
  else
  if (storage_scheme.compare(MM_SKEW_STR) == 0)
    mm_set_skew(&(ccm->matcode));
  else
    return MM_UNSUPPORTED_TYPE;      

  // check matcode
  if ((ccm->matcode)[3] != 'S') {
    printf("Only symmetric matrices are supported.\n");
    exit(1);
  }

  // read past comment block
  do {
    getline(input_file,line);
  } while (line[0] == '%');

  // get matrix dimensions
  iss.clear();
  iss.str(line);
  iss >> (ccm->m) >> (ccm->n) >> (ccm->nnz);

  // initialize specint
  ccm->a = -1.0;
  ccm->b = 1.0;

  // allocate host memory
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
    getline(input_file,line);
    iss.clear();
    iss.str(line);
    if ((ccm->matcode)[2] != 'P') {
      iss >> (ccm->rowinds)[ii] >> (ccm->colinds)[ii] >> (ccm->vals)[ii];
    }
    else {
      iss >> (ccm->rowinds)[ii] >> (ccm->colinds)[ii];
      (ccm->vals)[ii] = 1.0;
    }
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
  input_file.close();

  // create cublas handle
  if(cublasCreate(&(ccm->cublashandle)) != 0) {
    printf("CUBLAS initialization failed.\n");
    exit(1);
  }

  // set pointer mode to Host
  cublasSetPointerMode(ccm->cublashandle,CUBLAS_POINTER_MODE_HOST);

  // create cusparse handle
  if(cusparseCreate(&(ccm->cusparsehandle)) != 0) {
    printf("CUSPARSE initialization failed.\n");
    exit(1);
  }

  // set pointer mode to Host
  cusparseSetPointerMode(ccm->cusparsehandle,CUSPARSE_POINTER_MODE_HOST);

  // create cusparse MatDescr
  if(cusparseCreateMatDescr(&(ccm->matdescr)) != 0) {
    printf("CUSPARSE MatDescr initialization failed.\n");
    exit(1);
  }

  // allocate device memory
  if(cudaMalloc(&(ccm->drowinds),((ccm->m)+1)*sizeof(int)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccm->dcolinds),(ccm->nnz)*sizeof(int)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccm->dvals),(ccm->nnz)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }
  if(cudaMalloc(&(ccm->dtemp),2*(ccm->m)*sizeof(double)) != 0) {
    printf("Memory allocation failed.\n");
    exit(1);
  }

  // sort entries of CCM
  cuchebmatrix_sort(ccm);

  // convert CCM to csr format and copy to GPU
  cuchebmatrix_csr(ccm);

  // return  
  return 0;

}

