#include <cucheb.h>

/* driver */
int main(){

  // local variables
  double e1, e2;
  double vecs[4];

  // call symeig
  cuchebutils_2x2symeig(4,4,1e-15,&e1,&e2,&vecs[0]);

  // print eigs
  printf( " e1 = %+e\n",e1);
  printf( " e2 = %+e\n",e2);

  // print vecs
  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<2; jj++) {
      printf( " %+e",vecs[2*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // return 
  return 0;

}
