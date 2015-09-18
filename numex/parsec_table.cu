#include <cucheb.h>

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("./numex/parsec_data1.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/parsec_data1.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{c|c|c|c|c|c|c}\n";
  output_file << "Matrix & $n$ & $nnz$";
  //output_file << " & $bsize$ & $nblocks$ & $niters$ & $ndotprods$";
  //output_file << " & $degree$";
  output_file << " & $matvecs$";
  //output_file << " & preproc. (sec) & Lanczos (sec)";
  output_file << " & time (s)";
  output_file << " & $nconv$";
  output_file << " & $maxres$ \\\\\\hline\n";

  // variables to parse file
  string line, matname;
  int n, nnz, bsize, nblocks, niters, ndotprods, maxdeg, nmatvecs, nconv;
  double preproc, lanczos, maxres;

  // scan first line
  getline(input_file,line);

  // loop through lines
  while (!input_file.eof()) {

    // read in data
    input_file >> matname >> n >> nnz >> bsize >> nblocks >> niters >> 
                  ndotprods >> maxdeg >> nmatvecs >> preproc >> lanczos >>
                  nconv >> maxres;

    // exit if end of file
    if(input_file.eof()) { break; }

    // write to file
    output_file << matname << " & ";
    output_file << "$" << n << "$" << " & ";
    output_file << "$" << nnz << "$" << " & ";
    //output_file << "$" << bsize << "$" << " & ";
    //output_file << "$" << nblocks << "$" << " & ";
    //output_file << "$" << niters << "$" << " & ";
    //output_file << "$" << ndotprods << "$" << " & ";
    //output_file << "$" << maxdeg << "$" << " & ";
    output_file << "$" << nmatvecs << "$" << " & ";
    //output_file << "$" << preproc << "$" << " & ";
    //output_file << "$" << lanczos << "$" << " & ";
    output_file << "$" << fixed << setprecision(0) << preproc+lanczos << "$" << " & ";
    output_file << "$" << nconv << "$" << " & ";
    output_file << "$" << scientific << setprecision(1) << maxres << "$" << " \\\\\n";

  }

  // output_file footer
  output_file << "\\end{tabular}";

  // close input file
  input_file.close();

  // close output file
  output_file.close();





  // attempt to open input file
  input_file.open("./numex/parsec_data2.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/parsec_data2.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{c|c|c|c|c|c|c}\n";
  output_file << "Matrix & $n$ & $nnz$";
  //output_file << " & $bsize$ & $nblocks$ & $niters$ & $ndotprods$";
  //output_file << " & $degree$";
  output_file << " & $matvecs$";
  //output_file << " & preproc. (sec) & Lanczos (sec)";
  output_file << " & time (s)";
  output_file << " & $nconv$";
  output_file << " & $maxres$ \\\\\\hline\n";

  // scan first line
  getline(input_file,line);

  // loop through lines
  while (!input_file.eof()) {

    // read in data
    input_file >> matname >> n >> nnz >> bsize >> nblocks >> niters >> 
                  ndotprods >> maxdeg >> nmatvecs >> preproc >> lanczos >>
                  nconv >> maxres;

    // exit if end of file
    if(input_file.eof()) { break; }

    // write to file
    output_file << matname << " & ";
    output_file << "$" << n << "$" << " & ";
    output_file << "$" << nnz << "$" << " & ";
    //output_file << "$" << bsize << "$" << " & ";
    //output_file << "$" << nblocks << "$" << " & ";
    //output_file << "$" << niters << "$" << " & ";
    //output_file << "$" << ndotprods << "$" << " & ";
    //output_file << "$" << maxdeg << "$" << " & ";
    output_file << "$" << nmatvecs << "$" << " & ";
    //output_file << "$" << preproc << "$" << " & ";
    //output_file << "$" << lanczos << "$" << " & ";
    output_file << "$" << fixed << setprecision(0) << preproc+lanczos << "$" << " & ";
    output_file << "$" << nconv << "$" << " & ";
    output_file << "$" << scientific << setprecision(1) << maxres << "$" << " \\\\\n";

  }

  // output_file footer
  output_file << "\\end{tabular}";

  // close input file
  input_file.close();

  // close output file
  output_file.close();
  // return
  return 0;

}
