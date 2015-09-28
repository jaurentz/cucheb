#include <cucheb.h>

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("./numex/parsec/parsec_data.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/tables/parsec_table.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{l|c|c|c|c|c|c}\n";
  output_file << "Matrix & degree & iters & matvecs & dotprods & time & residual ";

  // variables to parse file
  string line, matname;
  int n, nnz, bsize, nblocks, niters, ndotprods, maxdeg, nmatvecs, nconv;
  double preproc, lanczos, maxres;

  int exponent;
  double mantissa;

  // scan first line
  getline(input_file,line);

  // loop through lines
  while (!input_file.eof()) {

    // read 3 at a time
    for (int ii=0; ii<3; ii++) {

      // read in data
      input_file >> matname >> n >> nnz >> bsize >> nblocks >> niters >> 
                    ndotprods >> maxdeg >> nmatvecs >> preproc >> lanczos >>
                    nconv >> maxres;

      // exit if end of file
      if(input_file.eof()) { break; }

      // write to file

      // newline line
      if (ii==0) { output_file << "\\\\\\hline\n"; }
      else { output_file << "\\\\\n"; }

      // matrix name
      if (ii==1) { output_file << "\\verb|" << matname << "|"; }
      else { output_file << ""; }

      // degree
      if (maxdeg > 99) { output_file << " & $" << maxdeg << "$"; }
      else { output_file << " & $\\phantom{0}" << maxdeg << "$"; }

      // niters
      if (niters > 9) { output_file << " & $" << niters << "$"; }
      else { output_file << " & $\\phantom{0}" << niters << "$"; }

      // nmatvecs
      if (nmatvecs > 99999) { output_file << " & $" << nmatvecs << "$"; }
      else if (nmatvecs > 9999) { output_file << " & $\\phantom{0}" << nmatvecs << "$"; }
      else if (nmatvecs > 999) { output_file << " & $\\phantom{00}" << nmatvecs << "$"; }
      else if (nmatvecs > 99) { output_file << " & $\\phantom{000}" << nmatvecs << "$"; }
      else { output_file << " & $\\phantom{0000}" << nmatvecs << "$"; }

      // ndotprods
      if (ndotprods > 99999) { output_file << " & $" << ndotprods << "$"; }
      else if (ndotprods > 9999) { output_file << " & $\\phantom{0}" <<
                                                  ndotprods << "$"; }
      else if (ndotprods > 999) { output_file << " & $\\phantom{00}" <<
                                                 ndotprods << "$"; }
      else if (ndotprods > 99) { output_file << " & $\\phantom{000}" <<
                                                ndotprods << "$"; }
      else { output_file << " & $\\phantom{0000}" << ndotprods << "$"; }

      // time
      if (preproc+lanczos >= 100 ) { output_file << " & $" << fixed <<
                                                     setprecision(0) << preproc+lanczos << "$"; }
      else if (preproc+lanczos >= 10 ) { output_file << " & $\\phantom{0}" << fixed <<
                                                         setprecision(0) << preproc+lanczos << "$"; }
      else { output_file << " & $\\phantom{00}" << fixed << setprecision(0) << preproc+lanczos << "$"; }

      // maxres
      exponent = floor(log10(abs(maxres)));
      mantissa = maxres/pow(10.0,exponent);
      output_file << " & $" << fixed << setprecision(1) << mantissa <<
                     "e{-" << abs(exponent) << "}$ ";

    }

  }

  // final line
  output_file << " \\\\\n"; 

  // output_file footer
  output_file << "\\end{tabular}";

  // close input file
  input_file.close();

  // close output file
  output_file.close();

  // return
  return 0;

}
