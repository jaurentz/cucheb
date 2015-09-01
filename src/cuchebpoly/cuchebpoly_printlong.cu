#include <cuchebpoly.h>

/* routine for long print */
int cuchebpoly_printlong(cuchebpoly* ccp){

  // print banner
  printf("\ncuchebpoly:\n");

  // degree
  printf(" degree = %d\n",ccp->degree);
 
  // a and b
  printf(" [a,b] = [%+e,%+e]\n",ccp->a,ccp->b);
  
  // coeffs
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    printf(" coeffs[%d] = %+e\n",ii,ccp->coeffs[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

