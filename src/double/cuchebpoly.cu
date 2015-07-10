#include <cuchebpoly.h>

/* routine for basic initialization */
int cuchebpoly_init(cuchebpoly* ccp){

  // set degree
  ccp->degree = 0;
  
  // set a and b
  ccp->a = -1.0;
  ccp->b = 1.0;
  
  // set coeffs
  ccp->coeffs[0] = 1.0;
  for (int ii=0; ii<DOUBLE_DEG; ii++) {
    ccp->coeffs[ii+1] = 0.0;
  }

  // return 
  return 0;

}

/* routine for standard print */
int cuchebpoly_print(cuchebpoly* ccp){

  // print banner
  printf("\ncuchebpoly:\n");

  // degree
  printf(" degree = %d\n",ccp->degree);
 
  // a and b
  printf(" [a,b] = [%+e,%+e]\n",ccp->a,ccp->b);
  
  // coeffs
  for (int ii=0; ii<ccp->degree+1; ii++) {
    printf(" coeffs[%d] = %+e\n",ii,ccp->coeffs[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

/* routine for long print */
int cuchebpoly_printlong(cuchebpoly* ccp){

  // print banner
  printf("\ncuchebpoly:\n");

  // degree
  printf(" degree = %d\n",ccp->degree);
 
  // a and b
  printf(" [a,b] = [%+e,%+e]\n",ccp->a,ccp->b);
  
  // coeffs
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    printf(" coeffs[%d] = %+e\n",ii,ccp->coeffs[ii]);
  }
  printf("\n");

  // return 
  return 0;

}

/* routine for Chebyshev points */
int cuchebpoints(double a, double b, double* points){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set points
  double alpha = (b-a)/2.0;
  double beta = (b+a)/2.0;
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    points[ii] = alpha*sin(DOUBLE_PI*(2.0*ii-DOUBLE_DEG)/(2.0*DOUBLE_DEG)) +
                  beta;
  }

  // return 
  return 0;

}

/* cuchebcoeffs */
int cuchebcoeffs (double *coeffs){
 
  // allocate workspace
  cufftDoubleReal *input;
  cudaMalloc(&input,2*DOUBLE_DEG*sizeof(cufftDoubleReal));
  cufftDoubleComplex *output;
  cudaMalloc(&output,(DOUBLE_DEG+1)*sizeof(cufftDoubleComplex));
 
  // initialize cufft
  cufftHandle cufftHand;
  cufftPlan1d(&cufftHand, 2*DOUBLE_DEG, CUFFT_D2Z, 1);
 
  // initialize input 
  int deg = DOUBLE_DEG;
  int N = 2*deg;

  cudaMemcpy(&input[0], &coeffs[deg], sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(&input[deg], &coeffs[0], sizeof(double), cudaMemcpyHostToDevice);
  for (int ii=1; ii<deg; ii++) {
    cudaMemcpy(&input[ii], &coeffs[deg-ii], sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(&input[N-ii], &coeffs[deg-ii], sizeof(double), cudaMemcpyHostToDevice);
  }

  // execute plan
  cufftExecD2Z(cufftHand,input,output);
 
  // extract output
  cudaMemcpy(&coeffs[0], &output[0], sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&coeffs[deg], &output[deg], sizeof(double), cudaMemcpyDeviceToHost);
  for (int ii=1; ii<deg; ii++) {
    cudaMemcpy(&coeffs[ii], &output[ii], sizeof(double), cudaMemcpyDeviceToHost);
  }

  // normalize output
  coeffs[0] = coeffs[0]/(double)(deg)/2.0;
  coeffs[deg] = coeffs[deg]/(double)(deg)/2.0;
  for (int ii=1; ii<deg; ii++) {
    coeffs[ii] = coeffs[ii]/(double)(deg);
  }

  // free cufft
  cufftDestroy(cufftHand);
 
  // free workspace
  cudaFree(input);
  cudaFree(output);
 
  // return success
  return 0;

}

/* routine for chopping Chebyshev coefficients */
int cuchebchop(int* degree, double* coeffs){

  // find maximum
  double maximum = 0.0;
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    maximum = max(maximum, abs(coeffs[ii]));
  }

  // maximum == 0
  if (maximum == 0) {
    *degree = 0;
  }
  
  // maximum != 0
  else {
    
    // compute degree
    for (int ii=0; ii<MAX_DOUBLE_DEG+1; ii++) {

      // set current degree
      *degree = MAX_DOUBLE_DEG-ii;

      // exit if trailing coefficient is too large
      if (abs(coeffs[MAX_DOUBLE_DEG-ii]) >= DOUBLE_EPS*maximum) {
        break;
      }

    }

  }

  // return 
  return 0;

}

/* routine for creating point filter */
int cuchebpoly_pointfilter(double a, double b, double rho, cuchebpoly* ccp){

  // check a and b
  if ( a >= b ) {
    return 1;
  }

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // compute Chebyshev points in [a,b]
  cuchebpoints(a,b,&(ccp->coeffs[0]));

  // compute shift
  double shift;
  if (rho <= a) {shift = a;}
  else if (rho >= b) {shift = b;}
  else {shift = rho;}

  // compute function values for f(x) = exp(-100*(x-shift)^2)
  double scl = pow(b - a,2);
  for (int ii=0; ii<DOUBLE_DEG+1; ii++) {
    ccp->coeffs[ii] = exp(-100.0*pow(ccp->coeffs[ii]-shift,2)/scl);
  }
 
  // compute Chebyshev coefficients
  cuchebcoeffs(&(ccp->coeffs[0]));

  // chop Chebyshev coefficients
  cuchebchop(&(ccp->degree),&(ccp->coeffs[0]));

  // return 
  return 0;

}

/* routine for creating step filter */
int cuchebpoly_stepfilter(double a, double b, double c, double d, cuchebpoly* ccp){

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
  ccp->degree = MAX_DOUBLE_DEG;

  // set a and b in ccp
  ccp->a = a;
  ccp->b = b;

  // scale everything to [-1,1]
  lb = (2.0*lb - (b+a))/(b-a);
  ub = (2.0*ub - (b+a))/(b-a);
 
  // compute Chebyshev coefficients
  double pi = DOUBLE_PI;
  double aclb = acos(lb);
  double acub = acos(ub);
  ccp->coeffs[0] = (aclb - acub)/pi;
  for (int ii=1; ii<MAX_DOUBLE_DEG+1; ii++) {
    ccp->coeffs[ii] = 2.0*(sin(ii*aclb) - sin(ii*acub))/(ii*pi);
  }

  // apply Jackson damping
  int deg = MAX_DOUBLE_DEG;
  double alpha = 1.0/(deg+2.0);
  double beta = sin(pi*alpha);
  double gamma = cos(pi*alpha);
  for (int ii=0; ii<MAX_DOUBLE_DEG+1; ii++) {
    ccp->coeffs[ii] = alpha*((deg+2.0-ii)*beta*cos(ii*pi*alpha) +
                       sin(ii*pi*alpha)*gamma)*ccp->coeffs[ii]/beta;
  }

  // chop Chebyshev coefficients
  cuchebchop(&(ccp->degree),&(ccp->coeffs[0]));

  // return 
  return 0;

}

