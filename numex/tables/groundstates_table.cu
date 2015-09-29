#include <cucheb.h>

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("./numex/groundstates/groundstates_data.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/tables/groundstates_table.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{l|c|c|c|c|c|c}\n";
  output_file << "\\hline\n";
  output_file << "\\multirow{2}{*}{Matrix}" <<
                 " & \\multirow{2}{*}{eigs} & \\multirow{2}{*}{deg}" <<
                 " & \\multirow{2}{*}{iters} & \\multirow{2}{*}{matvecs}" <<
                 " & \\multirow{2}{*}{time}" <<
                 " & \\multirow{2}{*}{residual} \\\\\n";
  output_file << " & & & & & & \\\\\\hline\n";
  output_file << "\\hline\n";

  // variables to parse file
  string matname;
  int n, nnz, bsize, nblocks, niters, ndotprods, maxdeg, nmatvecs, nconv;
  double a, b, neigs, preproc, lanczos, maxres;

  int exponent;
  double mantissa;
  int ones, thousands;

  // loop through lines
  while (!input_file.eof()) {

    // read 3 at a time
    for (int ii=0; ii<3; ii++) {

      // read in data
      input_file >> matname >> a >> b >> neigs >> n >> nnz >> bsize >> nblocks >> 
                    niters >> ndotprods >> maxdeg >> nmatvecs >> preproc >> lanczos >>
                    nconv >> maxres;

      // exit if end of file
      if(input_file.eof()) { break; }

      // write to file

      // matrix name
      if (ii==1) { output_file << "\\verb|" << matname << "|"; 

        // neigs
        if (nconv > 99) { output_file << " & $" << setprecision(0) << nconv << "$"; }
        else { output_file << " & $\\phantom{0}" << setprecision(0) << nconv << "$"; }
      }
      else { output_file << " & &"; }

      // degree
      if (maxdeg > 99) { output_file << " & $" << setprecision(0) << maxdeg << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << maxdeg << "$"; }

      // niters
      if (niters > 9) { output_file << " & $" << setprecision(0) << niters << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << niters << "$"; }

      // nmatvecs
      ones = nmatvecs%1000;
      thousands = (nmatvecs-ones)/1000;
      if (thousands > 99) {
        output_file << " & $" << thousands << 
                       "," << setw(3) << setfill('0') << setprecision(3) << ones << "$";
      }
      else {
        output_file << " & $\\phantom{0}" << thousands << 
                       "," << setw(3) << setfill('0') << setprecision(3) << ones << "$";
      }

      // ndotprods
//      ones = ndotprods%1000;
//      thousands = (ndotprods-ones)/1000;
//      if (thousands > 99) {
//        output_file << " & $" << thousands << 
//                       "," << setw(3) << setfill('0') << setprecision(3) << ones << "$";
//      }
//      else {
//        output_file << " & $\\phantom{0}" << thousands << 
//                       "," << setw(3) << setfill('0') << setprecision(3) << ones << "$";
//      }

      // time
      if (preproc+lanczos >= 100 ) { output_file << " & $" << fixed <<
                                                     setprecision(0) << preproc+lanczos << "$"; }
      else if (preproc+lanczos >= 10 ) { output_file << " & $\\phantom{0}" << fixed <<
                                                         setprecision(0) << preproc+lanczos << "$"; }
      else { output_file << " & $\\phantom{00}" << fixed << setprecision(0) << 
                            preproc+lanczos << "$"; }

      // maxres
      exponent = floor(log10(abs(maxres)));
      mantissa = maxres/pow(10.0,exponent);
      output_file << " & $" << fixed << setprecision(1) << mantissa <<
                     "e{-" << abs(exponent) << "}$ ";

      // new line
      if (ii==2) { output_file << "\\\\\\hline\n"; }
      else { output_file << "\\\\\n"; }

    }

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
