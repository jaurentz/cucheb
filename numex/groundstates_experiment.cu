#include <cucheb.h>

/* driver */
int main(){

  // set device
  //cudaSetDevice(1);

  // matrix files root directory
  const string rootdir("../matrices/");

  // number of matrices
  const int nummats = 7;

  // number of trials
  const int numtrials = 3;

  // matrix names
  string matnames[nummats*numtrials] = { "Ge87H76",
                                         "Ge87H76",
                                         "Ge87H76",
                                         "Ge99H100",
                                         "Ge99H100",
                                         "Ge99H100",
                                         "Si41Ge41H72",
                                         "Si41Ge41H72",
                                         "Si41Ge41H72",
                                         "Qdot3",
                                         "Qdot3",
                                         "Qdot3",
                                         "Si87H76",
                                         "Si87H76",
                                         "Si87H76",
                                         "Ga41As41H72",
                                         "Ga41As41H72",
                                         "Ga41As41H72",
                                         "Laplacian",
                                         "Laplacian",
                                         "Laplacian" };

  // matrix names
  int degrees[nummats][numtrials] = { {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {50,100,-1},
                                      {100,200,-1},
                                      {200,300,-1} };


  // output file
  string ofile("./numex/groundstates_data.txt" );

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
      cuchebmatrix_expertlanczos(100, -10.0, degrees[ii][jj],
                                 3, DEF_NUM_VECS, DEF_STEP_SIZE,
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
