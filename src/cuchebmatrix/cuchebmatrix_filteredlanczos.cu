#include <cucheb.h>
/*
  cuchebmatrix_init

  This routine computes the all the eigenvalues of a matrix in user-specified
  interval using the filtered Lanczos procedure. The following inputs are
  required:

    lbnd  - the lower bound of the desired interval
    ubnd  - the upper bound of the desired interval
    bsize - the size of the Lanczos blocks
    ccm   - a reference to an initialized instance of a cuchebmatrix
    ccl   - a reference to an uninitialized instance of a cucheblanczos

*/

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
