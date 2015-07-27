/*-------------------*/
/*   Matrix Formats  */
/*-------------------*/

/*------- COO format */
typedef struct coo_t_ {
  int n;
  int nnz;
  int *ir;
  int *jc;
  double *val;
} coo_t;

/*------- CSR format */
typedef struct csr_t_ {
  int n;
  int nnz;
  int *ia;
  int *ja;
  double *a;
} csr_t;

/*------- JAD format */
typedef struct jad_t_ {
  int n;
  int nnz;
  int *ia;
  int *ja;
  double *a;
  int njad;
  int *perm;
  double *w;
} jad_t;

/*------- DIA format */
typedef struct dia_t_ {
  int n;
  int nnz;
  int ndiags;
  int stride;
  double *diags;
  int *ioff;
} dia_t;


/*--------------------------*
        Matrix Wrapper
 *--------------------------*/

typedef struct matrix_t_ {

/*------------------------- *
  n is the size of matrix
  nnz is the number of non-zeros
  NOTE: nnz may be changed in 
  csr/jad->nnz for padding zeros
 *------------------------- */

  int n, nnz;
/*-------- csr */
  csr_t *h_csr, *d_csr;
/*-------- jad */
  jad_t *h_jad, *d_jad;
/*-------- dia */
  dia_t *h_dia, *d_dia;
/*---- spmv */
  void (*spmv) (struct matrix_t_*, double*, double*, int);
} matrix_t;

