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

  // scale to [-1,1]
  double A, B;
  A = max((2.0*lb - (b+a))/(b-a),-1.0);
  B = min((2.0*ub - (b+a))/(b-a),1.0);

  // local variables
  int degree;
  //degree = ceil(3.0*(DOUBLE_PI)/abs(acos(A)-acos(B)));
  degree = ceil(8.0*(DOUBLE_PI)/abs(acos(A)-acos(B)));

  // create stepfilter
  cuchebpoly_stepfilter(a,b,lb,ub,degree,ccp);

  // return 
  return 0;

}

