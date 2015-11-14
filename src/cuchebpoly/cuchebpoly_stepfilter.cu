#include <cucheb.h>

/* routine for creating step filter */
int cuchebpoly_stepfilter(double a, double b, double c, double d, int order, cuchebpoly* ccp){

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

  // save transition points
  double p1, p2;
  p1 = lb;
  p2 = ub;

  // scale everything to [-1,1]
  lb = max((2.0*lb - (b+a))/(b-a),-1.0);
  ub = min((2.0*ub - (b+a))/(b-a),1.0);
 
  // compute Chebyshev coefficients
  int deg;
  deg = ccp->degree;
  double pi = DOUBLE_PI;
  double aclb = acos(lb);
  double acub = acos(ub);
  ccp->coeffs[0] = (aclb - acub)/pi;
  for (int ii=1; ii<deg+1; ii++) {
    ccp->coeffs[ii] = 2.0*(sin(ii*aclb) - sin(ii*acub))/(ii*pi);
  }

  // apply Jackson damping
//  double alpha = 1.0/(deg+2.0);
//  double beta = sin(pi*alpha);
//  double gamma = cos(pi*alpha);
//  for (int ii=0; ii<deg+1; ii++) {
//    ccp->coeffs[ii] = alpha*((deg+2.0-ii)*beta*cos(ii*pi*alpha) +
//                       sin(ii*pi*alpha)*gamma)*ccp->coeffs[ii]/beta;
//  }

  // compute transition values
  p1 = cuchebpoly_clenshaw(ccp,p1); 
  p2 = cuchebpoly_clenshaw(ccp,p2); 

  // compute scale factor
  p1 = .5/min(abs(p1),abs(p2));

  // scale coefficients
  for (int ii=0; ii<deg+1; ii++) {
    ccp->coeffs[ii] = p1*ccp->coeffs[ii];
  }

  // return 
  return 0;

}

