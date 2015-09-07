#include <cucheb.h>

/* driver */
int main(){

  // local variables
  double e1, e2;
  double c, s;

  // call symeig
  cuchebutils_2x2symeig(4,4,1e-15,&e1,&e2,&c,&s);

  // print eigs
  printf( " e1 = %+e\n",e1);
  printf( " e2 = %+e\n",e2);

  // print vecs
  printf( " %+e, %+e\n",c,-s);
  printf( " %+e, %+e\n",s,c);
  printf("\n");

  // return 
  return 0;

}
