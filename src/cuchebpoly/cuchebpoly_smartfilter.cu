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

  // set degree
  ccp->degree = 0;

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // scale everything to [-1,1]
  lb = max((2.0*lb - (b+a))/(b-a),-1.0);
  ub = min((2.0*ub - (b+a))/(b-a),1.0);
 
  // set some local variables for computing coefficients
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

  // initialize norm
  double tol, nrm, nrm_exact;
  nrm = pi*(ccp->coeffs[0]*ccp->coeffs[0]);
  nrm_exact = sqrt(aclb-acub);
  //tol = pow(2.0,-4)*nrm_exact;
  tol = 5.0e-2*nrm_exact;
  //tol = pow(2.0,-5)*nrm_exact;
  
  // reduce tolerance if subinterval contains an end point
  if (lb == -1.0 || ub == 1.0 ){ tol = .25*tol; }

  // loop through degrees
  for (int jj=0; jj<MAX_DOUBLE_DEG; jj++) {

    // update coefficients
    ccp->coeffs[jj+1] = 2.0*(sin((jj+1)*aclb) - sin((jj+1)*acub))/((jj+1)*pi);
    deg += 1;
    ccp->degree = deg;
    nrm += pi*(ccp->coeffs[jj+1]*ccp->coeffs[jj+1])/2.0;

    // update transition points
    p1 += ccp->coeffs[jj+1]*cos((double)(jj+1)*aclb);
    p2 += ccp->coeffs[jj+1]*cos((double)(jj+1)*acub);

    // check error
    //if (abs(sqrt(nrm)-nrm_exact) < tol && jj >= 9) { break; }
    if (abs(sqrt(nrm)-nrm_exact) < tol) { break; }

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

