#include <cucheb.h>

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  string temp;
  string rootdir("/home/aurentz/Projects/CUCHEB/cucheb/numex/");
  string matdir("/home/aurentz/Projects/CUCHEB/matrices/");
  ifstream input_file;
  ofstream output_file;
  cuchebmatrix ccm;

  // attempt to open output file
  temp = rootdir + "tables/matrices_info.txt";
  output_file.open( temp.c_str() );
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // variables to parse file
  string matname;

  // attempt to open input file
  temp = rootdir + "tables/matrices.txt";
  input_file.open( temp.c_str() );
  if (!input_file.is_open()) { 
    printf("Could not open matrix file.\n");
    exit(1); 
  }

  // loop through lines
  while (!input_file.eof()) {

    // read in data
    input_file >> matname;

    // exit if end of file
    if(input_file.eof()) { break; }

    // initialize matrix
    temp = matdir + matname + ".mtx";
    cuchebmatrix_init(temp, &ccm);

    // compute spectral interval
    cuchebmatrix_specint(&ccm);

    // write to file
    output_file << matname << " " << ccm.m << " " << ccm.nnz << " " << 
                   setprecision(15) << ccm.a << " " << setprecision(15) << ccm.b << "\n";

    // destroy cuchebmatrix
    cuchebmatrix_destroy(&ccm);

  }

  // close input file
  input_file.close();

  // close output file
  output_file.close();

  // return
  return 0;

}
