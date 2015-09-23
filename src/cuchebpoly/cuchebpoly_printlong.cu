#include <cucheb.h>

/* routine for long print */
int cuchebpoly_printlong(cuchebpoly* ccp){

  // print banner
  printf("\ncuchebpoly:\n");

  // degree
  printf(" degree = %d\n",ccp->degree);
 
  // a and b
  printf(" [a,b] = [%+e,%+e]\n",ccp->a,ccp->b);
  printf("\n");
  
  // points
//  for (int ii=0; ii < 2*MAX_DOUBLE_DEG; ii++) {
//    printf(" points[%d] = %+e\n",ii,ccp->points[ii]);
//  }
//  printf("\n");
  
  // coeffs
  for (int ii=0; ii<ccp->degree+1; ii++) {
    printf(" coeffs[%d] = %+e\n",ii,ccp->coeffs[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

