#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // matrix files root directory
  const string rootdir("../matrices/DIMACS/");

  // number of matrices
  const int nummats = 10;

  // number of trials
  const int numtrials = 3;

  // matrix names
  string matnames[nummats*numtrials] = { "144",
                                         "144",
                                         "144",
                                         "auto",
                                         "auto",
                                         "auto",
                                         "ca2010",
                                         "ca2010",
                                         "ca2010",
                                         "caidaRouterLevel",
                                         "caidaRouterLevel",
                                         "caidaRouterLevel",
                                         "coPapersDBLP",
                                         "coPapersDBLP",
                                         "coPapersDBLP",
                                         "delaunay_n20",
                                         "delaunay_n20",
                                         "delaunay_n20",
                                         "fe_ocean",
                                         "fe_ocean",
                                         "fe_ocean",
                                         "m14b",
                                         "m14b",
                                         "m14b",
                                         "mn2010",
                                         "mn2010",
                                         "mn2010",
                                         "rgg_n_2_20_s0",
                                         "rgg_n_2_20_s0",
                                         "rgg_n_2_20_s0" };

  // filter degrees
  int degrees[nummats][numtrials] = { {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1} };

  // block sizes
  int bsize[nummats] = { 1,
                         2,
                         1,
                         1,
                         1,
                         1,
                         1,
                         1,
                         1,
                         2 };


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
      cuchebmatrix_expertlanczos(50, 1.0e300, degrees[ii][jj],
                                 bsize[ii], DEF_NUM_VECS, DEF_STEP_SIZE,
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
