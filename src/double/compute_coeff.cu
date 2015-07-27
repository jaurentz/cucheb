
#include "gpusollib.h"

void compute_coeff(int m, double a, double b, int damping, double* mu)
{

   // <math.h> is included in "gpusollib.h"

   /* Set up initial variables */

   double alpha1 = a;
   double alpha2 = b;

   double thetJ = PI / (m+2); // PI is defined in "gpusollib.h"
   double thetL = PI / (m+1);

   double a1 = 1 / (m+2);
   double a2 = sin(thetJ);

   double beta1 = acos(alpha1);
   double beta2 = acos(alpha2);

   /* Main loop to compute the polynomial coefficients */

   double jac;
   int k;

   for (k = 0; k <= m; k++) {
     
       // Depends on what the initial damping is
       if (damping == 0) {
          jac = 1;
       } else if (damping == 1) {
	 jac = (a1*sin((k+1)*thetJ)/a2) + (1-(k+1)*a1)*cos(k*thetJ);
       } else if (damping == 2) {
          jac = 1;
          if (k > 0)
             jac = sin(k*thetL) / (k*thetL);          
       }

       if (k == 0)
          mu[k] = -jac*( beta2-beta1 ) / PI;
       else
          mu[k] = -2*jac*( sin(k*beta2)-sin(k*beta1) ) / (PI*k);

   }

   // Print coefficients -- for debugging purposes only
   //for ( k = 0; k <= m; k++) 
   //   printf("%d'th coefficient: %f\n", k, mu[k]);

   /* Done computing the coefficients - now check validity 

   double *vkm1   = (double*) malloc(n*sizeof(double));
   double *vk     = (double*) malloc(n*sizeof(double));
   double *vkp1   = (double*) malloc(n*sizeof(double));
   double *result = (double*) malloc(n*sizeof(double));
   double *yi     = (double*) malloc(n*sizeof(double));

   memset(vk, 1.0, n*sizeof(double));

   int scal;

   for (int k = 0; k <= m; k++) {

       int one = 1;
       serial_daxpy(n, mu[k], vk, yi);
       scal = 2;

       if (k==0)
         scal = 1;
       
       vector_operations(3, n, xi, vk, result);
       serial_dscal(n, scal, result);
       vector_operations(2, n, result, vkm1, vkp1);
       memcpy(vk, vkm1, n*sizeof(double));
       memcpy(vkp1, vk, n*sizeof(double));
   }
   
   */

}













