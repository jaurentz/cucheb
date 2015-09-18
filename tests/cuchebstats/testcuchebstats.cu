#include <cucheb.h>

/* driver */
int main(){

  // output file
  string ofile("/home/aurentz/Desktop/data.csv");

  // number of matrices
  const int nummats = 3;

  // matname string array
  string matnames[nummats] = { "Tom", "Dick", "Harry" };

  // cuhebstats
  cuchebstats ccstats[nummats];

  // print to file
  cuchebstats_fileprint(ofile,nummats,&matnames[0],&ccstats[0]);

  // return 
  return 0;

}
