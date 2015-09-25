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
  if (neig > DEF_NUM_VECS) {
    printf("\ncuchebmatrix_filteredlanczos:\n");
    printf(" Number of desired eigenvalues is too large!\n\n");
    exit(1);
  }

  // make sure shift is valid
  if (isnan(shift)) {
    printf("\ncuchebmatrix_filteredlanczos:\n");
    printf(" Shift cannot be NaN!\n\n");
    exit(1);
  }

  // call expert lanczos 
  cuchebmatrix_expertlanczos(neig,shift,-1,
                             bsize,DEF_NUM_VECS,DEF_STEP_SIZE,
                             ccm,ccl,ccstats);

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

  // call expert lanczos
  cuchebmatrix_expertlanczos(lbnd, ubnd, -1,
                             bsize, DEF_NUM_VECS, DEF_STEP_SIZE,
                             ccm, ccl, ccstats);

  // return  
  return 0;

}
