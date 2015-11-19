#include <cucheb.h>

/* routine for creating step filter */
int cuchebpoly_stepfilter(double a, double b, double c, double d, int order, 
                          cuchebpoly* ccp){

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

  // set degree
  ccp->degree = min(max(0,order),MAX_DOUBLE_DEG);

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // scale everything to [-1,1]
  lb = max((2.0*lb - (b+a))/(b-a),-1.0);
  ub = min((2.0*ub - (b+a))/(b-a),1.0);
 
  // compute Chebyshev coefficients
  int deg;
  deg = ccp->degree;
  double pi = DOUBLE_PI;
  double aclb = acos(lb);
  double acub = acos(ub);

  // initialize coefficients
  ccp->coeffs[0] = (aclb - acub)/pi;

  // save transition points
  double p1, p2;
  p1 = ccp->coeffs[0];
  p2 = p1;

  // loop through degrees
  for (int jj=0; jj<deg; jj++) {

    // update coefficients
    ccp->coeffs[jj+1] = 2.0*(sin((jj+1)*aclb) - sin((jj+1)*acub))/((jj+1)*pi);

    // update transition points
    p1 += ccp->coeffs[jj+1]*cos((double)(jj+1)*aclb);
    p2 += ccp->coeffs[jj+1]*cos((double)(jj+1)*acub);

  }

  // compute scale factor
  if (lb == -1.0) { p1 = .5/abs(p2); }
  else if (ub == 1.0) { p1 = .5/abs(p1); }
  else { p1 = .5/min(abs(p1),abs(p2)); }

  // scale coefficients
  for (int ii=0; ii<deg+1; ii++) {
    ccp->coeffs[ii] = p1*ccp->coeffs[ii];
  }

  // return 
  return 0;

}

