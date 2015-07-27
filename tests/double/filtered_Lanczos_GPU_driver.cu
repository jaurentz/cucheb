#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <datatype.h>
#include <protos.h>

/*-------------------------------*/
/*        Filtered Lanczos       */
/*-------------------------------*/

int main() {

   //----------
   // Variables
   //----------
   char       finput1[] = "matrix.mtx";
   char       finput2[] = "input.txt";
   FILE*      input_h;
   matrix_t   *mat;
   coo_t      *h_coo;
   int        n, nnz, k;

   int    lanczos_steps, degree, mat_choice;
   double lmin, lmax;

   //------------------
   //GPU initialization
   //------------------
   /*
   printf("%d\n", gpusol_init);
   if (gpusol_init()) {
      printf("GPU-Solv Init Error\n");
      return 1;
   }
   */

   //----------------------------------------
   // Read the matrix file and allocate space
   //----------------------------------------
   // The matrix is stored in COO (Matrix-Market) format
   Calloc(mat, 1, matrix_t);
   Calloc(h_coo, 1, coo_t);
   read_coo_MM(h_coo, finput1);
   n   = h_coo->n;   // number of rows
   nnz = h_coo->nnz; // number of non-zero entries

   //----------------------------------------
   // Read input file to determine parameters
   //----------------------------------------
   input_h = fopen(finput2, "r");
   for (k = 0; k < 5; k++) {
      if (k==0)
        fscanf(input_h, "%d", &mat_choice);
      else if (k==1)
        fscanf(input_h, "%d", &lanczos_steps);
      else if (k==2)
        fscanf(input_h, "%d", &degree);
      else if (k==3)
        fscanf(input_h, "%lg", &lmin);
      else if (k==4)
        fscanf(input_h, "%lg", &lmax);
   } 

   for (k = 0; k < 5; k++) {
      if (k==0)
        printf("%d\n", mat_choice);
      else if (k==1)
        printf("%d\n", lanczos_steps);
      else if (k==2)
        printf("%d\n", degree);
      else if (k==3)
        printf("%f\n", lmin);
      else if (k==4)
        printf("%f\n", lmax);
   } 


   // Initialize all possible matrix-structures
   mat->n   = n;  
   mat->nnz = nnz;
   Calloc( mat->h_csr, 1, csr_t );
   Calloc( mat->h_jad, 1, jad_t );
   Calloc( mat->h_dia, 1, dia_t );
   
   // CSR matrix convertion
   if (mat_choice == 0)
      COO2CSR( h_coo, mat->h_csr );

   setup_matrix( mat, mat_choice );

   //----------------------------------
   // Lanczos with polynomial filtering
   //----------------------------------
   lanczos_steps  = 100;      // maximum number of Lanczos steps
   degree  = 200;             // degree of the polynomial
   lmin = 0;                  // for now lambda_min should be known
   lmax = 13;                 // same applies to lambda_max
   double w    = (lmax-lmin) / 2; // scaling to [-1,1]
   double c    = (lmax+lmin) / 2; //        >>

   filtered_lanczos(mat, lanczos_steps, degree, w, c, mat_choice);

   //--------------------------------
   // Done, free the allocated memory
   //--------------------------------
   //free_matrix(mat);
   //free_coo(h_coo);

   //-------------------------
   // Finalize and check errors
   //-------------------------
   //gpusol_finalize();

}





