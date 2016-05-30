#include <cucheb.h>

/* driver */
int main(){

  // compute variables
  string mtx("");
  cuchebmatrix ccm;
  cucheblanczos ccl;
  cuchebstats ccstats;

  // initialize matrix
  cuchebmatrix_init(mtx, &ccm);

  // call filtered lanczos for an interval
  cuchebmatrix_expertlanczos(a, b, deg, bsize, nvecs, ssize,
                               &ccm, &ccl, &ccstats);

  // print stats
  cuchebstats_print(&ccstats);

  // destroy CCL
  cucheblanczos_destroy(&ccl);

  // destroy cuchebmatrix
  cuchebmatrix_destroy(&ccm);

  // return 
  return 0;

}
