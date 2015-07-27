
#include "gpusollib.h"

// Routine that performs entry-wise addition/subtraction/multiplication of two vectors
// No size checking


void vector_operations(int choice, int n, double *x, double *y, double *z) {

  if (choice == 1) { // addition

   for (int i = 0; i < n; i++)
      z[i] = x[i] + y[i];

  }

  else if (choice == 2) { // subtraction

   for (int i = 0; i < n; i++)
      z[i] = x[i] - y[i];
  
  }
  else if (choice == 3) { // multiplication

   for (int i = 0; i < n; i++)
       z[i] = x[i] * y[i];

  }
}

void serial_daxpy(int n, double a, double *x, double *y) {


   for (int i = 0; i < n; i++)
       y[i] = y[i] + a*x[i];


}


void serial_dscal(int n, double a, double *x) {

   for (int i =0; i < n; i++)
     x[i] = a*x[i];

}
