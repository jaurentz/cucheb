# include <stdio.h>
# include <cuda.h>
# include "magma.h"
# include "magma_lapack.h"

int main (int argc , char ** argv){

	magma_init (); // initialize Magma
	magma_int_t n=1024 , n2=n*n;
	double *a, *r; // a, r - nxn matrices on the host
	double *d_r ; // nxn matrix on the device
	double * h_work ; // workspace
	magma_int_t lwork ; // h_work size
	magma_int_t * iwork ; // workspace
	magma_int_t liwork ; // iwork size
	double *w1 , *w2; // w1 ,w2 - vectors of eigenvalues
	double error , work [1]; // used in difference computations
	magma_int_t ione = 1 , i, j, info ;
	double mione = -1.0;
	magma_int_t incr = 1 , inci = 1;
	magma_dmalloc_cpu (&w1 ,n); // host memory for real
	magma_dmalloc_cpu (&w2 ,n); // eigenvalues
	magma_dmalloc_cpu (&a,n2 ); // host memory for a
	magma_dmalloc_cpu (&r,n2 ); // host memory for r
	magma_dmalloc (& d_r ,n2 ); // device memory for d_r
	
	// Query for workspace sizes
	double aux_work[1];
	magma_int_t aux_iwork[1];
	
	magma_dsyevd_gpu('V','L',n,d_r,n,w1,r,n,aux_work,-1,aux_iwork,-1,&info);
	lwork =(magma_int_t)aux_work[0];
	liwork = aux_iwork[0];
	iwork =(magma_int_t*)malloc(liwork*sizeof(magma_int_t));
	magma_dmalloc_cpu(&h_work,lwork); // memory for workspace
	
	// define a, r // [1 0 0 0 0 ...]
	for(i=0;i<n;i ++){ // [0 2 0 0 0 ...]
		a[i*n+i]=( double )(i +1); // a = [0 0 3 0 0 ...]
		r[i*n+i]=( double )(i +1); // [0 0 0 4 0 ...]
	} // [0 0 0 0 5 ...]
	printf (" upper left corner of a:\n"); // .............
	magma_dprint (5 ,5 ,a,n); // print part of a
	magma_dsetmatrix ( n, n, a, n, d_r , n); // copy a -> d_r
	
	// compute the eigenvalues and eigenvectors for a symmetric,
	// real nxn matrix ; Magma version
	magma_dsyevd_gpu(MagmaVec,MagmaLower,n,d_r,n,w1,r,n,h_work,lwork,iwork,liwork,&info);
	printf (" first 5 eigenvalues of a:\n");
	for(j=0;j <5;j++) printf ("%f\n",w1[j]); // print first eigenvalues
	printf (" left upper corner of the matrix of eigenvectors :\n");
	magma_dgetmatrix ( n, n, d_r , n, r, n ); // copy d_r -> r
	magma_dprint (5 ,5 ,r,n); // part of the matrix of eigenvectors
	
	free (w1 ); // free host memory
	free (w2 ); // free host memory
	free (a); // free host memory
	free (r); // free host memory
	free ( h_work ); // free host memory
	magma_free (d_r ); // free device memory
	magma_finalize (); // finalize Magma
	
	return EXIT_SUCCESS ;
}

// upper left corner of a:
//[
// 1.0000 0. 0. 0. 0.
// 0. 2.0000 0. 0. 0.
// 0. 0. 3.0000 0. 0.
// 0. 0. 0. 4.0000 0.
// 0. 0. 0. 0. 5.0000
// ];
// first 5 eigenvalues of a:
// 1.000000
// 2.000000
// 3.000000
// 4.000000
// 5.000000
// left upper corner of the matrix of eigenvectors :
//[
// 1.0000 0. 0. 0. 0.
// 0. 1.0000 0. 0. 0.
// 0. 0. 1.0000 0. 0.
// 0. 0. 0. 1.0000 0.
// 0. 0. 0. 0. 1.0000
// ];
// difference in eigenvalues : 0.000000 e+00



