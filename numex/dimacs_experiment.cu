#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // matrix files root directory
  const string rootdir("../matrices/DIMACS/");

  // number of matrices
  const int nummats = 14;

  // number of trials
  const int numtrials = 2;

  // matrix names
  string matnames[nummats*numtrials] = { "144",
                                         "144",
                                         "598a",
                                         "598a",
                                         "al2010",
                                         "al2010",
                                         "ar2010",
                                         "ar2010",
                                         "auto",
                                         "auto",
                                         "az2010",
                                         "az2010",
                                         "ca2010",
                                         "ca2010",
                                         "caidaRouterLevel",
                                         "caidaRouterLevel",
                                         "citationCiteseer",
                                         "citationCiteseer",
                                         "co2010",
                                         "co2010",
                                         "coAuthorsCiteseer",
                                         "coAuthorsCiteseer",
                                         "coAuthorsDBLP",
                                         "coAuthorsDBLP",
                                         "coPapersCiteseer",
                                         "coPapersCiteseer",
                                         "coPapersDBLP",
                                         "coPapersDBLP" };

  // matrix names
  int degrees[nummats][numtrials] = { {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100},
                                      {50,100} };


  // output file
  string ofile("./numex/dimacs_data.txt" );

  // cuchebstats array
  cuchebstats ccstats[nummats*numtrials]; 
  
  // local variables
  string mtxfile;
  cuchebmatrix ccm;
  cucheblanczos ccl;

  // loop through matrices
  for (int ii=0; ii<nummats; ii++) {

    // print matname
    printf(" %s",matnames[ii*numtrials].c_str());

    // set mtxfile
    mtxfile = rootdir + matnames[ii*numtrials] + ".mtx";

    // initialize matrix
    cuchebmatrix_init(mtxfile, &ccm);

    // trials with various degrees
    for (int jj=0; jj<numtrials; jj++) {

      // call filtered lanczos for an interval
      cuchebmatrix_expertlanczos(100, 1.0e300, degrees[ii][jj],
                                 1, DEF_NUM_VECS, DEF_STEP_SIZE,
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
