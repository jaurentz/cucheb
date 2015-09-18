#include <cucheb.h>

/* driver */
int main(){

  // matrix files root directory
  const string rootdir("../matrices/");

  // number of matrices
  const int nummats = 21;

  // matrix names
  string matnames[nummats] = { "Si2",
                               "SiH4",
                               "SiNa",
                               "Na5",
                               "benzene",
                               "Si10H16",
                               "Si5H12",
                               "SiO",
                               "Ga3As3H12",
                               "GaAsH6",
                               "H2O",
                               "Si34H36",
                               "Ge99H100",
                               "Ge87H76",
                               "Ga10As10H30",
                               "Ga19As19H42",
                               "SiO2",
                               "Si41Ge41H72",
                               "CO",
                               "Si87H76",
                               "Ga41As41H72" };


  // number experiments
  const int numexs = 2;

  // output file
  string ofile[numexs] = { "./numex/parsec_data1.txt",
                           "./numex/parsec_data2.txt" };

  // cuchebstats array
  cuchebstats ccstats[nummats]; 
  
  // local variables
  string mtxfile;
  cuchebmatrix ccm;
  cucheblanczos ccl;


  // Experiment 1
  // loop through matrices
  for (int ii=0; ii<nummats; ii++) {

    // print matname
    printf(" %s",matnames[ii].c_str());

    // set mtxfile
    mtxfile = rootdir + matnames[ii] + ".mtx";

    // initialize matrix
    cuchebmatrix_init(mtxfile, &ccm);

    // call filtered lanczos for a point
    cuchebmatrix_filteredlanczos(10, -1e100, 3, &ccm, &ccl, &ccstats[ii]);

    // print stats
    cuchebstats_print(&ccstats[ii]);

    // destroy CCM
    cuchebmatrix_destroy(&ccm);

    // destroy CCL
    cucheblanczos_destroy(&ccl);

  }

  // print ccstats to file
  cuchebstats_fileprint(ofile[0],nummats,&matnames[0],&ccstats[0]);




  // Experiment 2
  // loop through matrices
  for (int ii=0; ii<nummats; ii++) {

    // print matname
    printf(" %s",matnames[ii].c_str());

    // set mtxfile
    mtxfile = rootdir + matnames[ii] + ".mtx";

    // initialize matrix
    cuchebmatrix_init(mtxfile, &ccm);

    // call filtered lanczos for an interval
    cuchebmatrix_filteredlanczos(-10.0, 0.0, 3, &ccm, &ccl, &ccstats[ii]);

    // print stats
    cuchebstats_print(&ccstats[ii]);

    // destroy CCM
    cuchebmatrix_destroy(&ccm);

    // destroy CCL
    cucheblanczos_destroy(&ccl);

  }

  // print ccstats to file
  cuchebstats_fileprint(ofile[1],nummats,&matnames[0],&ccstats[0]);

  // return 
  return 0;

}
