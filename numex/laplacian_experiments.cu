#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // matrix files root directory
  const string rootdir("../matrices/");

  // number of matrices
  const int nummats = 1;

  // lower bounds
  double lb[nummats] = { 1.0 };

  // upper bounds
  double ub[nummats] = { 1.01 };

  // number of trials
  const int numtrials = 3;

  // matrix names
  string matnames[nummats*numtrials] = { "Laplacian",
                                         "Laplacian",
                                         "Laplacian" };

  // matrix names
  int degrees[nummats][numtrials] = { {1000,1600,-1} };


  // output file
  string ofile("./numex/laplacian_data.txt" );

  // cuchebstats array
  cuchebstats ccstats[nummats*numtrials]; 
  
  // local variables
  string mtxfile;
  cuchebmatrix ccm;
  cucheblanczos ccl;

  // loop through matrices
  for (int ii=0; ii<nummats; ii++) {

    // print matname
    printf(" %s",matnames[ii].c_str());

    // set mtxfile
    mtxfile = rootdir + matnames[ii] + ".mtx";

    // initialize matrix
    cuchebmatrix_init(mtxfile, &ccm);

    // trials with various degrees
    for (int jj=0; jj<numtrials; jj++) {

      // call filtered lanczos for an interval
      cuchebmatrix_expertlanczos(lb[ii], ub[ii], degrees[ii][jj],
                                 3, 1440, DEF_STEP_SIZE,
                                 &ccm, &ccl, &ccstats[ii*numtrials+jj]);

      // print stats
      cuchebstats_print(&ccstats[ii*numtrials+jj]);

      // destroy CCL
      cucheblanczos_destroy(&ccl);

    }

    // destroy CCM
    cuchebmatrix_destroy(&ccm);

  }

  // print ccstats to file
  cuchebstats_fileprint(ofile,nummats*numtrials,&matnames[0],&ccstats[0]);

  // return 
  return 0;

}
