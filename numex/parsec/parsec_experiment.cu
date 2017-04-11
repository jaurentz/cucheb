#include <cucheb.h>

/* driver */
int main(){

  // set device
  cudaSetDevice(1);

  // compute variables
  string temp;
//  string rootdir("/LUSTRE/users/jaurentz/Projects/CUCHEB/cucheb/numex/");
//  string matdir("/LUSTRE/users/jaurentz/Projects/CUCHEB/matrices/");
  string rootdir("/home/saady/kalantzi/CUCHEB/cucheb/numex/");
  string matdir("/home/saady/kalantzi/CUCHEB/matrices/");
  ifstream input_file;
  ofstream output_file;
  cuchebmatrix ccm;
  cucheblanczos ccl;
  cuchebstats ccstats;

  // attempt to open output file
  temp = rootdir + "parsec/parsec_data.txt";
  output_file.open( temp.c_str() );
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // variables to parse file
  string matname;
  double a, b;
  int neigs, deg, bsize, nvecs, ssize;

  // attempt to open input file
  temp = rootdir + "parsec/parsec_matrices.txt";
  input_file.open( temp.c_str() );
  if (!input_file.is_open()) { 
    printf("Could not open matrix file.\n");
    exit(1); 
  }

  // loop through lines
  while (!input_file.eof()) {

    // loop through sub experiments
    for (int jj=0; jj<5; jj++){ 

      // read in data
      input_file >> matname >> a >> b >> neigs >> deg >> bsize >> nvecs >> ssize;
  
      // exit if end of file
      if(input_file.eof()) { break; }
  
      // print matrix name
      printf("\nMatrix: %s\n",matname.c_str());
  
      // initialize matrix
      temp = matdir + matname + ".mtx";
      if (jj==0){ cuchebmatrix_init(temp, &ccm); }
  
      // call filtered lanczos for an interval
      cuchebmatrix_expertlanczos(a, b, deg, bsize, nvecs, ssize,
                                   &ccm, &ccl, &ccstats);
  
      // print stats
      cuchebstats_print(&ccstats);
  
      // write to file
      output_file << matname.c_str() << " "; 
      output_file << setprecision(15) << a << " ";
      output_file << setprecision(15) << b << " ";
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
      output_file << ccstats.innerprod_time << " ";
      output_file << ccstats.matvec_time << " ";
      output_file << ccstats.total_time << " ";
      output_file << ccstats.num_conv << " ";
      output_file << ccstats.max_res << "\n";
  
      // destroy CCL
      cucheblanczos_destroy(&ccl);
  
      // destroy cuchebmatrix
      if (jj==4) { cuchebmatrix_destroy(&ccm); }

    }

  }

  // close input file
  input_file.close();

  // close output file
  output_file.close();

  // return 
  return 0;

}
