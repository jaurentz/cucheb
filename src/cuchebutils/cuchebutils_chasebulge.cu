#include <cucheb.h>

/* function to perform banded symmetric bulge chase */
int cuchebutils_chasebulge(int n, int bwidth, double* bands, int ldbands,
                           double* bulge, double* vecs, int ldvecs){

  // local variables
  int nsteps, ind;
  double c, s, nrm, a1, a2, b;

  // set nsteps
  nsteps = ceil((double)(n-bwidth)/(bwidth-1));
//printf(" nsteps = %d\n\n", nsteps);

  // loop for chasing
  for (int ii=0; ii < nsteps; ii++) {

    // set ind
    ind = ii*ldbands*(bwidth-1) + bwidth - 1;

    // generate rotation
    cuchebutils_rotation(bands[ind],*bulge,&c,&s,&nrm);
  
    // remove bulge
    bands[ind] = c*bands[ind] + s*(*bulge);

    // update rows
    for (int jj=0; jj < bwidth-2; jj++) {

      // update ind
      ind = ind + ldbands - 1;

      // update bands
      b = c*bands[ind] + s*bands[ind+1];
      bands[ind+1] = -s*bands[ind] + c*bands[ind+1];
      bands[ind] = b;

    }

    // update core 2x2
    ind = ind + ldbands - 1;
    a1 = bands[ind];
    b = bands[ind+1];
    a2 = bands[ind+ldbands];

    bands[ind] = c*(a1*c + b*s) + s*(b*c + a2*s);
    bands[ind+1] = c*(b*c - a1*s) + s*(a2*c - b*s);
    bands[ind+ldbands] = c*(a2*c - b*s) - s*(b*c - a1*s);

    // update cols
    ind = ind + ldbands + 1;
    for (int jj=0; jj < bwidth-2; jj++) {

      // update bands
      b = c*bands[ind-ldbands+1+jj] + s*bands[ind+jj];
      bands[ind+jj] = -s*bands[ind-ldbands+1+jj] + c*bands[ind+jj];
      bands[ind-ldbands+1+jj] = b;

    }

    // update bulge
    ind = ind + bwidth-2;
    *bulge = s*bands[ind];
    bands[ind] = c*bands[ind];

    // update transforming matrix
    ind = (ii+1)*ldvecs*(bwidth-1);
    for (int jj=0; jj < ldvecs; jj++) {
 
      b = c*vecs[ind+jj] + s*vecs[ind+ldvecs+jj];
//printf(" b = %+e\n", b);
      vecs[ind+ldvecs+jj] = -s*vecs[ind+jj] + c*vecs[ind+ldvecs+jj];
      vecs[ind+jj] = b;

    }
  
  } 

  // return success
  return 0;

}
