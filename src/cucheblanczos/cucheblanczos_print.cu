#include <cucheb.h>

/* routine for standard print */
int cucheblanczos_print(cucheblanczos* ccl){

  // print banner
  printf("\ncucheblanczos:\n");

  // print n
  printf(" n = %d\n",ccl->n);
 
  // print bsize
  printf(" bsize = %d\n",ccl->bsize);
 
  // print nblocks
  printf(" nblocks = %d\n",ccl->nblocks);
  printf("\n");
 
  // return 
  return 0;

}

