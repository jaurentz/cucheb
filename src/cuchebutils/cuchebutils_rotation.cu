#include <cucheb.h>

/* function to generate givens rotation */
int cuchebutils_rotation(const double a, const double b, double* c, double* s,
                         double* nrm){

  // initialize c, s and nrm
  *c = 1.0;
  *s = 0.0;
  *nrm = 0.0;

  // return if a and/or b is NAN
  if ( a != a || b != b ){
    printf("\ncuchebutils_rotation:\n");
    printf(" a and b must not be NAN!\n\n");
    exit(1);
  }

  // compute rotation
  if ( abs(a) == 0.0 && abs(b) == 0.0 ) { return 0; }
  else if ( abs(a) >= abs(b) ){

    *s = b/a;
    *nrm = copysign(1.0,a)*sqrt(1.0 + (*s)*(*s));    
    *c = 1.0/(*nrm);
    *s = (*s)*(*c);
    *nrm = a*(*nrm);

  }
  else {

    *c = a/b;
    *nrm = copysign(1.0,b)*sqrt(1.0 + (*c)*(*c));    
    *s = 1.0/(*nrm);
    *c = (*s)*(*c);
    *nrm = b*(*nrm);

  }

  // return success
  return 0;

}
