#include <cucheb.h>

/* routine for standard print */
int cuchebblocklanczos_print(cuchebblocklanczos* ccb){

  // print banner
  printf("\ncuchebblocklanczos:\n");

  // print n
  printf(" n = %d\n",ccb->n);
 
  // print bsize
  printf(" bsize = %d\n",ccb->bsize);
 
  // print nblocks
  printf(" nblocks = %d\n",ccb->nblocks);
  printf("\n");
 
  // return 
  return 0;

}

