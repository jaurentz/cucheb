#include <cucheb.h>

/* y = p(x) */
double cuchebpoly_clenshaw(cuchebpoly* ccp, double x){

  // local variables
  int deg;
  double a, b;
  double* coeffs;
  deg = ccp->degree;
  a = ccp->a;
  b = ccp->b;
  coeffs = &(ccp->coeffs)[0];

  // scalars
  double y, dt1, dt2;
  double A, B;
  A = 4.0/(b-a);
  B = -2.0*(b+a)/(b-a);

  // initialize y
  y = coeffs[deg];
 
  // initialize dt1
  dt1 = 0.0;
 
  // loop for clenshaw
  for(int ii=0; ii<deg; ii++){

    // copy dt1 to dt2
    dt2 = dt1;

    // copy y to dt1
    dt1 = y;

    // scale A and B if ii == deg-1
    if(ii == deg-1){
      A = A/2.0;
      B = B/2.0;
    }

    // apply matrix
    y = A*x*dt1 + B*y;

    // add x
    y = coeffs[deg-ii-1] + y;

    // subtract dt2
    y = y - dt2;

  }

  // return 
  return y;

}

