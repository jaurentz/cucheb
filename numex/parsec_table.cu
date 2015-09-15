#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;


/* routine to convert csv file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("/Users/jared/Projects/CUCHEB/cucheb/numex/parsec_data.csv");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("/Users/jared/Projects/CUCHEB/cucheb/numex/parsec_data.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{ccccccccccccc}\n";
  output_file << "Matrix & $n$ & $nnz$ & $bsize$ & $nblocks$ & $niters$";
  output_file << " & $ndotprods$ & $maxdeg$ & $nmatvecs$ & preproc. (sec)";
  output_file << " & Lanczos (sec) & $nconv$ & $maxres$ \\\\\n";

  // variables to parse banner
//  string line, banner, mtx, crd, data_type, storage_scheme;

  // scan first line
//  getline(input_file,line);
//  istringstream iss(line);
//  iss >> banner >> mtx >> crd >> data_type >> storage_scheme;

  // output_file footer
  output_file << "\\end{tabular}";

  // close input file
  input_file.close();

  // close output file
  output_file.close();

  // return
  return 0;

}
