#include <cucheb.h>

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("./numex/dimacs/dimacs_data.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/tables/dimacs_table.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{l|c|c|c|c|c|c}\n";
  output_file << "\\hline\n";
  output_file << "\\multirow{2}{*}{Matrix}" << " & \\multirow{2}{*}{deg}" <<
                 " & \\multirow{2}{*}{iters} & \\multirow{2}{*}{matvecs}" <<
                 " & \\multirow{2}{*}{time} & \\multirow{2}{*}{eigs}" <<
                 " & \\multirow{2}{*}{residual} \\\\\n";
  output_file << " & & & & & & \\\\\\hline\n";
  output_file << "\\hline\n";

  // variables to parse file
  string matname;
  int n, nnz, bsize, nblocks, niters, ndotprods, maxdeg, nmatvecs, nconv;
  double a, b, preproc, innerprod, matvec, total, maxres;

  int exponent;
  double mantissa;
  int ones, thousands;

  // loop through lines
  while (!input_file.eof()) {

    // read in data
    input_file >> matname >> a >> b >> n >> nnz >> bsize >> nblocks >> 
                  niters >> ndotprods >> maxdeg >> nmatvecs >> preproc >> 
                  innerprod >> matvec >> total >> nconv >> maxres;

    // exit if end of file
    if(input_file.eof()) { break; }

    // write to file

    // matrix name
    output_file << "\\verb|" << matname << "|";

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

    // time
    if (total >= 100) { output_file << " & $" << fixed << setprecision(0) 
                                    << total << "$"; }
    else if (total >= 10) { output_file << " & $\\phantom{0}" << fixed 
                                        << setprecision(0) << total << "$"; }
    else { output_file << " & $\\phantom{00}" << fixed << setprecision(0) << 
                          total << "$"; }

    // nconv
    if (nconv > 99) { output_file << " & $" << setprecision(0) << nconv << "$"; }
    else if (nconv > 9) { output_file << " & $\\phantom{0}" << 
                                         setprecision(0) << nconv << "$"; }
    else { output_file << " & $\\phantom{00}" << setprecision(0) << nconv << "$"; }

    // maxres
    exponent = floor(log10(abs(maxres)));
    mantissa = maxres/pow(10.0,exponent);
    output_file << " & $" << fixed << setprecision(1) << mantissa;
    if (abs(exponent) > 9) { output_file << "e{-" << abs(exponent) << "}$ "; }
    else { output_file << "e{-0" << abs(exponent) << "}$ "; }

    // new line
    output_file << "\\\\\\hline\n";

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
