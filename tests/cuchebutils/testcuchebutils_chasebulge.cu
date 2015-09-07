#include <cucheb.h>

/* driver */
int main(){

  // create 4x4 banded symmetric matrix
  double bands[8] = { 2.0, -1.0, 
                      2.0, -1.0, 
                      2.0, -1.0, 
                      2.0,  0.0};

  // create matrix of vectors
  double vecs[16] = { 1.0, 0.0, 0.0, 0.0,
                      0.0, 1.0, 0.0, 0.0,
                      0.0, 0.0, 1.0, 0.0,
                      0.0, 0.0, 0.0, 1.0 };

  // print bands
  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<4; jj++) {
      printf( " %+e",bands[2*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // print vecs
  for (int ii=0; ii<4; ii++) {
    for (int jj=0; jj<4; jj++) {
      printf( " %+e",vecs[4*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // create bulge
  double bulge = 0.1;

  // chase bulge
  cuchebutils_chasebulge(4,2,&bands[0],2,&bulge,&vecs[0],4);

  // print bands
  for (int ii=0; ii<2; ii++) {
    for (int jj=0; jj<4; jj++) {
      printf( " %+e",bands[2*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // print vecs
  for (int ii=0; ii<4; ii++) {
    for (int jj=0; jj<4; jj++) {
      printf( " %+e",vecs[4*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // return 
  return 0;

}
