#include <cucheb.h>

/* reduce banded symmetric matrix to tridiagonal */
int cuchebutils_bandsymred(int n, int bwidth, double* bands, int ldbands,
                           double* vecs, int ldvecs){

  // local variables
  int ind;
  double a1, a2, b, c, s, bulge;

  // set local variables

  // outer loop to remove a band
  for (int ii=0; ii < bwidth-2; ii++) {
  //for (int ii=0; ii < 1; ii++) {

    // inner loop to remove subdiagonals
    for (int jj=0; jj < n-bwidth+1+ii; jj++) {
    //for (int jj=0; jj < 1; jj++) {

      // initialize bulge
      // create eliminator
      ind = jj*ldbands+bwidth-1-ii;
      cuchebutils_rotation(bands[ind-1],bands[ind],&c,&s,&bulge);
//printf( " c = %+e\n",c);
//printf( " s = %+e\n\n",s);

      // apply to bands
      // update rows
      ind = ind - 1;
      for (int kk=0; kk < bwidth-2-ii; kk++) {

        // update bands
        b = c*bands[ind] + s*bands[ind+1];
        bands[ind+1] = -s*bands[ind] + c*bands[ind+1];
        bands[ind] = b;

        // update ind
        ind = ind + ldbands - 1;

      }

      // update core
      a1 = bands[ind];
      b = bands[ind+1];
      a2 = bands[ind+ldbands];

      bands[ind] = c*(a1*c + b*s) + s*(b*c + a2*s);
      bands[ind+1] = c*(b*c - a1*s) + s*(a2*c - b*s);
      bands[ind+ldbands] = c*(a2*c - b*s) - s*(b*c - a1*s);
    
      // update cols
      ind = ind + ldbands + 1;
      for (int kk=0; kk < bwidth-2-ii; kk++) {

        // update bands
        b = c*bands[ind-ldbands+1+kk] + s*bands[ind+kk];
        bands[ind+kk] = -s*bands[ind-ldbands+1+kk] + c*bands[ind+kk];
        bands[ind-ldbands+1+kk] = b;

      }

      // compute bulge
      ind = ind+bwidth-2-ii;
      bulge = s*bands[ind];
      bands[ind] = c*bands[ind];
//printf( " bulge = %+e\n\n",bulge);

      // apply to vecs
      ind = (bwidth-2-ii+jj)*ldvecs;
      for (int kk=0; kk < ldvecs; kk++) {
 
        b = c*vecs[ind+kk] + s*vecs[ind+ldvecs+kk];
        vecs[ind+ldvecs+kk] = -s*vecs[ind+kk] + c*vecs[ind+ldvecs+kk];
        vecs[ind+kk] = b;

      }

      // chase bulge
      ind = bwidth-2-ii+jj;
      cuchebutils_chasebulge(n-ind,bwidth-ii,&bands[ind*ldbands],ldbands,
                             &bulge,&vecs[ind*ldvecs],ldvecs);

    }  

  } 

  // return success
  return 0;

}
