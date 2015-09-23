#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // matrix files root directory
  const string rootdir("../matrices/");

  // number of matrices
  const int nummats = 7;

  // matrix names
  string matnames[nummats] = { "Ge99H100",
                               "Ge87H76",
                               "Si41Ge41H72",
                               "Si87H76",
                               "Ga41As41H72",
                               "Andrews",
                               "Laplacian" };

  // lower bounds
  double lb[nummats] = { -0.645,
                         -0.650,
                         -0.640,
                         -0.660,
                         -0.640,
                          4.000,
                          1.000 };

  // upper bounds
  double ub[nummats] = { -0.00530,
                         -0.00960,
                         -0.00282,
                         -0.33000,
                         -0.00000,
                          5.00000,
                          1.01000 };


  // output file
  string ofile("./numex/parsec_data3.txt" );

  // cuchebstats array
  cuchebstats ccstats[nummats]; 
  
  // local variables
  string mtxfile;
  cuchebmatrix ccm;
  cucheblanczos ccl;


  // Experiment 3
  // loop through matrices
  for (int ii=0; ii<nummats; ii++) {

if (ii > 4) {break;}

    // print matname
    printf(" %s",matnames[ii].c_str());

    // set mtxfile
    mtxfile = rootdir + matnames[ii] + ".mtx";

    // initialize matrix
    cuchebmatrix_init(mtxfile, &ccm);
    cuchebmatrix_print(&ccm);

    // call filtered lanczos for an interval
    cuchebmatrix_filteredlanczos(lb[ii], ub[ii], 3, &ccm, &ccl, &ccstats[ii]);

    // print stats
    cuchebstats_print(&ccstats[ii]);

    // destroy CCM
    cuchebmatrix_destroy(&ccm);

    // destroy CCL
    cucheblanczos_destroy(&ccl);

  }

  // print ccstats to file
  cuchebstats_fileprint(ofile,nummats,&matnames[0],&ccstats[0]);

  // return 
  return 0;

}
