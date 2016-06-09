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
  output_file.open("./numex/tables/parsec_table2.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{l|c|c|c|c|c}\n";
  output_file << "\\hline\n";
  output_file << "\\multirow{2}{*}{Matrix} & \\multirow{2}{*}{$m$}" <<
                 " & \\multirow{2}{*}{iters} & \\multirow{2}{*}{PREPROC}" <<
                 " & \\multirow{2}{*}{ORTH}" <<
                 " & \\multirow{2}{*}{MV} \\\\\n";
  output_file << " & & & & & \\\\\\hline\n";
  output_file << "\\hline\n";

  // variables to parse file
  string matname;
  int neigs, n, nnz, bsize, nblocks, niters, ndotprods, maxdeg, nmatvecs, nconv;
  double a, b, preproc, innerprod, matvec, total, maxres;

  int ones;

  // loop through lines
  while (!input_file.eof()) {

    // read 5 at a time
    for (int ii=0; ii<5; ii++) {

      // read in data
      input_file >> matname >> a >> b >> neigs >> n >> nnz >> bsize >> nblocks >> 
                    niters >> ndotprods >> maxdeg >> nmatvecs >> preproc >> 
                    innerprod >> matvec >> total >> nconv >> maxres;

      // exit if end of file
      if(input_file.eof()) { break; }

      // write to file

      // matrix name
      if (ii==2) { output_file << "\\verb|" << matname << "|";}
      else { output_file << " "; }

      // degree
      if (maxdeg > 99) { output_file << " & $" << setprecision(0) << maxdeg << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << maxdeg << "$"; }

      // niters
      if (nblocks > 99) { output_file << " & $" << setprecision(0) << nblocks << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << nblocks << "$"; }

      // specint time
      ones = round(100.0*preproc/total);
      if (ones > 9) { output_file << " & $" << setprecision(0) << ones << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << ones << "$"; }
      output_file << "\\%";

      // innerprod time
      ones = round(100.0*innerprod/total);
      if (ones > 9) { output_file << " & $" << setprecision(0) << ones << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << ones << "$"; }
      output_file << "\\%";

      // matvec time
      ones = round(100.0*matvec/total);
      if (ones > 9) { output_file << " & $" << setprecision(0) << ones << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << ones << "$"; }
      output_file << "\\%";

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
