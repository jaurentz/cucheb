#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // compute variables
  string temp;
  string rootdir("/home/aurentz/Projects/CUCHEB/cucheb/numex/");
  string matdir("/home/aurentz/Projects/CUCHEB/matrices/");
  ifstream input_file;
  ofstream output_file;
  cuchebmatrix ccm;
  cucheblanczos ccl;
  cuchebstats ccstats;

  // attempt to open output file
  temp = rootdir + "groundstates/groundstates_data.txt";
  output_file.open( temp.c_str() );
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // variables to parse file
  string matname;
  double shift;
  int neigs, deg, bsize, nvecs, ssize;

  // attempt to open input file
  temp = rootdir + "groundstates/groundstates_matrices.txt";
  input_file.open( temp.c_str() );
  if (!input_file.is_open()) { 
    printf("Could not open matrix file.\n");
    exit(1); 
  }

  // loop through lines
  while (!input_file.eof()) {

    // read in data
    input_file >> matname >> shift >> neigs >> deg >> bsize >> nvecs >> ssize;

    // exit if end of file
    if(input_file.eof()) { break; }

    // initialize matrix
    temp = matdir + matname + ".mtx";
    cuchebmatrix_init(temp, &ccm);

    // call filtered lanczos for an interval
    cuchebmatrix_expertlanczos(neigs, shift, deg, bsize, nvecs, ssize,
                                 &ccm, &ccl, &ccstats);

    // print stats
    cuchebstats_print(&ccstats);

    // write to file
    output_file << matname.c_str() << " "; 
    output_file << neigs << " ";
    output_file << ccstats.mat_dim << " ";
    output_file << ccstats.mat_nnz << " ";
    output_file << ccstats.block_size << " ";
    output_file << ccstats.num_blocks << " ";
    output_file << ccstats.num_iters << " ";
    output_file << ccstats.num_innerprods << " ";
    output_file << ccstats.max_degree << " ";
    output_file << ccstats.num_matvecs << " ";
    output_file << ccstats.specint_time << " ";
    output_file << ccstats.arnoldi_time << " ";
    output_file << ccstats.num_conv << " ";
    output_file << ccstats.max_res << "\n";

    // destroy cuchebmatrix
    cuchebmatrix_destroy(&ccm);

    // destroy CCL
    cucheblanczos_destroy(&ccl);

  }

  // close input file
  input_file.close();

  // close output file
  output_file.close();

  // return 
  return 0;

}
