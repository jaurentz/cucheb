#include <cucheb.h>
/*
  cucheblanczos_print

  This routine prints some basic properties of an instance of a cucheblanczos
  object. The following inputs are required:

    ccl - a reference to an instance of a cucheblanczos

*/

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
 
  // print stop
  printf(" stop = %d\n",ccl->stop);
 
  // print nconv
  printf(" nconv = %d\n",ccl->nconv);
  printf("\n");
 
  // return 
  return 0;

}

