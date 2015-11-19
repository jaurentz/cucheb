#include <cucheb.h>


// replace substring in place
void ReplaceStringInPlace(std::string& subject, const std::string& search,
                          const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("./numex/engineering/engineering_data.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/tables/engineering_table2.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{l|c|c|c|c}\n";
  output_file << "\\hline\n";
  output_file << "\\multirow{2}{*}{Matrix} & \\multirow{2}{*}{deg}" <<
                 " & \\multirow{2}{*}{iters} & \\multirow{2}{*}{dotprods}" <<
                 " & \\multirow{2}{*}{matvecs} \\\\\n";
  output_file << " & & & & \\\\\\hline\n";
  output_file << "\\hline\n";

  // variables to parse file
  string matname;
  int n, nnz, bsize, nblocks, niters, ndotprods, maxdeg, nmatvecs, nconv;
  double a, b, per, preproc, innerprod, matvec, total, maxres;

  int ones;

  // loop through lines
  while (!input_file.eof()) {

    // read 3 at a time
    for (int ii=0; ii<2; ii++) {

      // read in data
      input_file >> matname >> a >> b >> per >> n >> nnz >> bsize >> nblocks >> 
                    niters >> ndotprods >> maxdeg >> nmatvecs >> preproc >> 
                    innerprod >> matvec >> total >> nconv >> maxres;

      // replace annoying underscores
      ReplaceStringInPlace(matname,"_","\\_");

      // exit if end of file
      if(input_file.eof()) { break; }

      // write to file

      // matrix name
      if (ii==0) { output_file << "\\multirow{2}{*}{\\texttt{" << matname << "}}"; }
      else { output_file << " "; }

      // degree
      if (maxdeg > 9) { output_file << " & $" << setprecision(0) << maxdeg << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << maxdeg << "$"; }

      // niters
      if (niters > 9) { output_file << " & $" << setprecision(0) << niters << "$"; }
      else { output_file << " & $\\phantom{0}" << setprecision(0) << niters << "$"; }
 
      // update total time
      if (ii==1) { total -= preproc; }

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
      if (ii==1) { output_file << "\\\\\\hline\n"; }
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
