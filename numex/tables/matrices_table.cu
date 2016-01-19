#include <cucheb.h>

/* routine to convert txt file to tex table */
int main(){

  // compute variables
  ifstream input_file;
  ofstream output_file;

  // attempt to open input file
  input_file.open("./numex/tables/matrices_info.txt");
  if (!input_file.is_open()) { 
    printf("Could not open input file.\n");
    exit(1); 
  }

  // attempt to open output file
  output_file.open("./numex/tables/matrices_table.tex");
  if (!output_file.is_open()) { 
    printf("Could not open output file.\n");
    exit(1); 
  }

  // output file banner
  output_file << "\\begin{tabular}{l|c|c|c|c}\n";
  output_file << "\\hline\n";
  output_file << "\\multirow{2}{*}{Matrix} & \\multirow{2}{*}{$n$}" <<
                 " & \\multirow{2}{*}{$nnz$} & \\multirow{2}{*}{$nnz/n$}" <<
                 " & \\multirow{2}{*}{Spectral interval} \\\\\n";
  output_file << " & & & & \\\\\\hline\n";
  output_file << "\\hline\n";

  // variables to parse file
  string matname;
  int n, nnz;
  double a, b;
  int ones, thousands, millions;
  int exponent;
  double mantissa;

  // loop through lines
  int ii = 0;
  while (!input_file.eof()) {

    // read in data
    input_file >> matname >> n >> nnz >> a >> b;

    // exit if end of file
    if(input_file.eof()) { break; }

    // write matname to file
    output_file << "\\verb|" << matname << "|";

    // write n to file
    ones = n%1000;
    thousands = ((n-ones)%1000000)/1000;
    millions = (n-ones-1000*thousands)/1000000;
    output_file << " & $";
    if (millions > 0) { output_file  << millions << ","; }
    else { output_file << "\\phantom{0,{}}"; }
    output_file << setw(3) << setfill('0') << setprecision(3) << thousands << "," << 
                   setw(3) << setfill('0') << setprecision(3) << ones << "$";

    // write nnz to file
    ones = nnz%1000;
    thousands = ((nnz-ones)%1000000)/1000;
    millions = (nnz-ones-1000*thousands)/1000000;
    output_file << " & $";
    if (millions > 9) { output_file  << millions << ","; }
    else if (millions > 0) { output_file  << "\\phantom{0}" <<  millions << ","; }
    else { output_file << "\\phantom{00,{}}"; }
    output_file << setw(3) << setfill('0') << setprecision(3) << thousands << "," << 
                   setw(3) << setfill('0') << setprecision(3) << ones << "$";

    // write nnz/n to file
    output_file << " & $";
    if (nnz/(1.0*n) > 9) { output_file << fixed << setprecision(1) << nnz/(1.0*n); }
    else { output_file  << "\\phantom{0}" << fixed << setprecision(1) 
                        << nnz/(1.0*n); }
    output_file << "$";

    // write specint to file
    output_file << " & $[";

    // write a
    exponent = floor(log10(abs(a)));
    mantissa = a/pow(10.0,exponent);
    if (mantissa >= 0) { output_file << "\\phantom{-{}}" << fixed <<
                         setprecision(2) << mantissa; }
    else { output_file << "-" << fixed << setprecision(2) << abs(mantissa); }
    if (exponent >= 0) { output_file << "e{+" << exponent << "},"; }
    else { output_file << "e{-" << abs(exponent) << "},"; }

    // write b
    exponent = floor(log10(abs(b)));
    mantissa = b/pow(10.0,exponent);
    if (mantissa >= 0) { output_file << "\\phantom{-{}}" << fixed <<
                         setprecision(2) << mantissa; }
    else { output_file << "-" << fixed << setprecision(2) << abs(mantissa); }
    if (exponent >= 0) { output_file << "e{+" << exponent << "}"; }
    else { output_file << "e{-" << abs(exponent) << "}"; }

    // new line
    if ( ii == 4 || ii == 8 ) { output_file << "]$ \\\\\\hline\n"; }
    else { output_file << "]$ \\\\\n"; }

    // update index
    ii++;

  }

  // output_file footer
    output_file << "\\hline\n";
  output_file << "\\end{tabular}";

  // close input file
  input_file.close();

  // close output file
  output_file.close();

  // return
  return 0;

}
