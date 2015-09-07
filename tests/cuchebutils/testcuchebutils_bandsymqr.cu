#include <cucheb.h>

/* driver */
int main(){

  // dimensions
  const int bwidth = 4;
  const int n = 6;
 
  // create banded symmetric matrix
  double bands[bwidth*n] = {0.0};
  for (int jj=0; jj<n; jj++) {
    bands[jj*bwidth] = 2.0;
  }
  for (int jj=0; jj<n-1; jj++) {
    bands[jj*bwidth+1] = 1.0;
  }
  for (int jj=0; jj<n-2; jj++) {
    bands[jj*bwidth+2] = 1.0;
  }
  for (int jj=0; jj<n-3; jj++) {
    bands[jj*bwidth+3] = 1.0;
  }

  // create matrix of vectors
  double vecs[n*n] = {0.0};
  for (int ii=0; ii<n; ii++) {
    for (int jj=0; jj<n; jj++) {
      if(ii==jj) {vecs[n*jj+ii] = 1.0;}
      else {vecs[jj*n+ii] = 0.0;}
    }
  }

  // print bands
  for (int ii=0; ii<bwidth; ii++) {
    for (int jj=0; jj<n; jj++) {
      printf( " %+e",bands[bwidth*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // print vecs
  for (int ii=0; ii<n; ii++) {
    for (int jj=0; jj<n; jj++) {
      printf( " %+e",vecs[n*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // array for evals
  double evals[n];

  // call bandsymqr
  cuchebutils_bandsymqr(n,bwidth,&bands[0],bwidth,&evals[0],&vecs[0],n);

  // print evals
  for (int ii=0; ii<n; ii++) {
    printf( " %+e",evals[ii]);
  }
  printf("\n\n");

  // print bands
  for (int ii=0; ii<bwidth; ii++) {
    for (int jj=0; jj<n; jj++) {
      printf( " %+e",bands[bwidth*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // print vecs
  for (int ii=0; ii<n; ii++) {
    for (int jj=0; jj<n; jj++) {
      printf( " %+e",vecs[n*jj+ii]);
    }
    printf("\n");
  }
  printf("\n");

  // return 
  return 0;

}
