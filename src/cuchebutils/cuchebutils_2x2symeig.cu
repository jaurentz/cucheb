#include <cucheb.h>

/* eigenvalues and eigenvectors of 2x2 symmetric matrix */
int cuchebutils_2x2symeig(double a1, double a2, double b, double* e1, double* e2,
                          double* c, double* s){

  // local variables
  double trace, detm, disc;
  trace = a1 + a2;
  detm = a1*a2 - b*b;
  disc = a1 - a2;
  disc = disc*disc + 4.0*b*b;
    
  // compute sqrt of discriminant 
  disc = sqrt(disc);
      
  // compute most accurate eigenvalue first
  if ( abs(trace+disc) > abs(trace-disc) ) {
    *e1 = (trace+disc)/2.0;
    *e2 = detm/(*e1);
  }
  else {
    *e1 = (trace-disc)/2.0;
    *e2 = detm/(*e1);
  }

  // compute eigenvectors
  cuchebutils_rotation(-b,a1-*e1,c,s,&trace);

  // return success
  return 0;

}
