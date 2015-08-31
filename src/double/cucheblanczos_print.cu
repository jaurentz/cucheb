#include <cucheblanczos.h>

/* routine for standard print */
int cucheblanczos_print(cucheblanczos* ccl){

  // print banner
  printf("\ncucheblanczos:\n");

  // print m
  printf(" n = %d\n",ccl->n);
 
  // print n
  printf(" nvecs = %d\n",ccl->nvecs);
  printf("\n");
 
  // return 
  return 0;

}

