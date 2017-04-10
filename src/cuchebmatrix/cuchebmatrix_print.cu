#include <cucheb.h>
/*
  cuchebmatrix_print

  This routine prints some basic properties of an instance of a cuchebmatrix
  object. The following inputs are required:

    ccm - a reference to an instance of a cuchebmatrix

*/

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
 
  // print [a,b]
  printf(" [a,b] = [%+e,%+e]\n",ccm->a,ccm->b);
  printf("\n");
 
  // return 
  return 0;

}

