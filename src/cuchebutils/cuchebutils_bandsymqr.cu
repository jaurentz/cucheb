#include <cucheb.h>

/* compute eigenvalues of banded symmetric matrix using qr */
int cuchebutils_bandsymqr(int n, int bwidth, double* bands, int ldbands,
                           double* evals, double* vecs, int ldvecs){

  // local variables
  int ind, start, stop, zero, maxits;
  double a1, a2, b, c, s, bulge, e1, e2;

  // set local variables
  start = 0;
  stop = n-1;
  zero = -1;
  maxits = 20*n;

  // reduce to tridiagonal
  cuchebutils_bandsymred(n, bwidth, bands, ldbands, vecs, ldvecs);

  // loop for chasing
  for (int ii=0; ii < maxits; ii++) {

    // exit if finished
    if ( stop < 0 ) { break; }

    // check for deflation
    for (int jj=0; jj < stop-start; jj++) {
      
      // tolerance
      ind = (stop-1-jj)*ldbands;
      b = DOUBLE_TOL*(abs(bands[ind])+abs(bands[ind+ldbands]));

      // check for zeros
      if (abs(bands[ind+1]) < b) {
        bands[ind+1] = 0.0;
        zero = stop-1-jj;
        start = zero+1;
        break;
      }

    }
//printf( " zero = %d\n",zero);
//printf( " start = %d\n",start);
//printf( " stop = %d\n\n",stop);

    // 1x1 case
    if (start == stop) {

      // store eval
      evals[stop] = bands[stop*ldbands];

      // update indices
      zero = -1;
      start = 0;
      stop = stop-1;

    }

    // 2x2 or greater
    else {

      // compute shift
      ind = (stop-1)*ldbands;
      a1 = bands[ind];
      b = bands[ind+1];
      a2 = bands[ind+ldbands];
//printf( " a1 = %+1.15e\n",a1);
//printf( " a2 = %+1.15e\n",a2);
//printf( " b = %+1.15e\n\n",b);


      cuchebutils_2x2symeig(a1,a2,b,&e1,&e2,&c,&s);
      bulge = bands[stop*ldbands];
      if ( abs(e1-bulge) < abs(e2-bulge) ){ bulge = e1; }
      else { bulge = e2; }
//printf( " e1 = %+1.15e\n",e1);
//printf( " e2 = %+1.15e\n\n",e2);

      // initialize bulge
      // create eliminator
      cuchebutils_rotation(bands[start*ldbands]-bulge,bands[start*ldbands+1],
                           &c,&s,&bulge);
//printf( " c = %+e\n",c);
//printf( " s = %+e\n\n",s);

      // apply to vecs
      ind = start*ldvecs;
      for (int jj=0; jj < ldvecs; jj++) {
 
        b = c*vecs[ind+jj] + s*vecs[ind+ldvecs+jj];
        vecs[ind+ldvecs+jj] = -s*vecs[ind+jj] + c*vecs[ind+ldvecs+jj];
        vecs[ind+jj] = b;

      }

      // apply to bands
      ind = start*ldbands;
      a1 = bands[ind];
      b = bands[ind+1];
      a2 = bands[ind+ldbands];

      bands[ind] = c*(a1*c + b*s) + s*(b*c + a2*s);
      bands[ind+1] = c*(b*c - a1*s) + s*(a2*c - b*s);
      bands[ind+ldbands] = c*(a2*c - b*s) - s*(b*c - a1*s);
    
      // compute bulge
      bulge = s*bands[ind+ldbands+1];
      bands[ind+ldbands+1] = c*bands[ind+ldbands+1];
//printf( " bulge = %+e\n\n",bulge);

      // chase bulge
      cuchebutils_chasebulge(stop-start+1,2,&bands[start*ldbands],ldbands,
                             &bulge,&vecs[start*ldvecs],ldvecs);

    }  

  // print bands
//  for (int kk=0; kk<2; kk++) {
//    for (int jj=0; jj<5; jj++) {
//      printf( " %+e",bands[2*jj+kk]);
//    }
//    printf("\n");
//  }
//  printf("\n");

  // print vecs
//  for (int kk=0; kk<5; kk++) {
//    for (int jj=0; jj<5; jj++) {
//      printf( " %+e",vecs[5*jj+kk]);
//    }
//    printf("\n");
//  }
//  printf("\n");

    // print error if not converged
    if (ii == maxits-1) {
      printf("cuchebutils_bandsymqr:\n");
      printf(" Did not converge!\n\n");
      return 1;
    }

  } 

  // return success
  return 0;

}
