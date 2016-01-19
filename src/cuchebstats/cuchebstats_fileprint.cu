#include <cucheb.h>

/* routine for printing to file */
int cuchebstats_fileprint(const string& fname, int nummats, string* matnames, cuchebstats* ccstats){

  // attempt to open file
  ofstream output_file; 
  output_file.open(fname.c_str());
  if (!output_file.is_open()) {
    printf("Could not open output file.\n");
    exit(1);
  }

  // print banner
  output_file << "mat_name mat_dim mat_nnz "; 
  output_file << "block_size num_blocks num_iters num_innerprods max_degree "; 
  output_file << "num_matvecs specint_time innerprod_time matvec_time total_time ";
  output_file << "num_conv max_res\n";

  // loop through data
  for (int ii=0; ii<nummats; ii++){

    // print matrix name
    output_file << matnames[ii].c_str() << " ";

    // print ccstats data
    output_file << ccstats[ii].mat_dim << " ";
    output_file << ccstats[ii].mat_nnz << " ";
    output_file << ccstats[ii].block_size << " ";
    output_file << ccstats[ii].num_blocks << " ";
    output_file << ccstats[ii].num_iters << " ";
    output_file << ccstats[ii].num_innerprods << " ";
    output_file << ccstats[ii].max_degree << " ";
    output_file << ccstats[ii].num_matvecs << " ";
    output_file << ccstats[ii].specint_time << " ";
    output_file << ccstats[ii].innerprod_time << " ";
    output_file << ccstats[ii].matvec_time << " ";
    output_file << ccstats[ii].total_time << " ";
    output_file << ccstats[ii].num_conv << " ";
    output_file << ccstats[ii].max_res << "\n";

  }

  // close file
  output_file.close();

  // return 
  return 0;

}

