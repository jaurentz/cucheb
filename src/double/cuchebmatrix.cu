#include <cuchebmatrix.h>

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
    getline(input_file,line);
    iss.clear();
    iss.str(line);
    iss >> (ccm->rowinds)[ii] >> (ccm->colinds)[ii] >> (ccm->vals)[ii];
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
  printf(" matcode = %.4s\n",ccm->matcode);
 
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
  printf(" matcode = %.4s\n",ccm->matcode);
 
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

/* routine for sorting entries */
int cuchebmatrix_sort(cuchebmatrix* ccm){

  // create a vector of pairs of pairs
  vector< pair< pair<int,int> , double > > mat;
  for(int ii=0; ii<(ccm->nnz); ii++){
    mat.push_back(make_pair( make_pair((ccm->rowinds)[ii],(ccm->colinds)[ii]) , (ccm->vals)[ii] ));
  }

  // sort vector
  sort(mat.begin(),mat.end());

  // update ccm
  for(int ii=0; ii < mat.size(); ii++){
    (ccm->rowinds)[ii] = mat[ii].first.first;
    (ccm->colinds)[ii] = mat[ii].first.second;
    (ccm->vals)[ii] = mat[ii].second;
  }

  // return 
  return 0;

}

/* routine for converting to csr format */
int cuchebmatrix_csr(cuchebmatrix* ccm){

  // loop through row inds
  int cind = 0;
  for(int ii=0; ii<(ccm->nnz); ii++){
    if((ccm->rowinds)[ii] > cind){
      cind += 1;
      (ccm->rowinds)[cind] = ii;
    }
  }
  (ccm->rowinds)[ccm->m] = ccm->nnz;

  // return 
  return 0;

}


