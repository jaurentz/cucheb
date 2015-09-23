#include <cucheb.h>

/* routine for creating smart filter */
int cuchebpoly_smartfilter(double a, double b, double c, double d, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // check c and d
  if ( c >= d ) {
    return 1;
  }

  // compute lower bound 
  double lb;
  if (c <= a) {lb = a;}
  else if (c >= b) {return 1;}
  else {lb = c;}

  // compute upper bound 
  double ub;
  if (d >= b) {ub = b;}
  else if (d <= a) {return 1;}
  else {ub = d;}

  // local variables
  int degree;
  degree = ceil(4.0*abs(b-a)/abs(ub-lb));
  //degree = min(degree,1000);

  // create stepfilter
  cuchebpoly_stepfilter(a,b,lb,ub,degree,ccp);

  // return 
  return 0;

}

