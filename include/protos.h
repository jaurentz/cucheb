int gpusol_init();

int gpusol_finalize();

double wall_timer();

void read_input( char* fn );

void COO2CSR(coo_t *coo, csr_t *csr);

void PadJAD32(jad_t *jad);

void CSR2JAD(csr_t *csr, jad_t *jad);

int CSR2DIA(csr_t *csr, dia_t *dia);

void csrcsc(int n, int n2, int job, int ipos, double *a, int *ja, int *ia, double *ao, int *jao, int *iao);

void sortrow(csr_t*);

void setup_matrix( matrix_t *mat, int choice );

int qsplitC(double *a, int *ind, int n, int ncut);

int read_coo_MM( coo_t *coo, char *filename );

void *cuda_malloc(int size);

void *cuda_malloc_host(int size);

void memcpyh2d(void *dest, void* src, int size);

void memcpyd2d(void *dest, void* src, int size);

void memcpyd2h(void *dest, void *src, int size);

void cuda_memset(void *, int, int);

void malloc_csr(int n, int nnz, csr_t *csr);

void realloc_csr(csr_t *csr, int nnz);

void cuda_malloc_csr(int n, int nnz, csr_t *d_csr);

void copy_csr_h2d(csr_t *h_csr, csr_t *d_csr);

void copy_csr_h2h(csr_t *csr1, csr_t *csr2);

void cuda_malloc_jad(int n, int njad, int nnz, jad_t *d_jad);

void copy_jad_h2d(jad_t *h_jad, jad_t *d_jad);

void cuda_malloc_dia(int nd, int strd, dia_t *d_dia);

void copy_dia_h2d(dia_t *h_dia, dia_t *d_dia);

void Free(void *p);

void cuda_free(void *p);

void cuda_free_host(void *p);

void free_coo(coo_t *coo);

void free_csr(csr_t *csr);

void cuda_free_csr(csr_t *d_csr);

void free_jad(jad_t *jad);

void cuda_free_jad(jad_t *d_jad);

void free_dia(dia_t *dia);

void cuda_free_dia(dia_t *d_dia);

void free_matrix(matrix_t *mat);

void scal_vec(int, double*, double *);

void spmv_csr_vector(matrix_t *mat, double *x, double *y, int neg);

void coo2csr(int n, int nnz, double *coo_vals, int *coo_rows, int *coo_cols, double *acsr, int *ia, int *ja);

//

void filtered_spmv_csr_vector(matrix_t *mat, double *x, double *y, int neg, int degree, double w, double c, double *mu, int mat_choice);

void vector_operations(int choice, int n, double *x, double *y, double *z);

void compute_coeff(int m, double a, double b, int damping, double* mu);

void serial_daxpy(int, double, double*, double*);

void serial_dscal(int, double, double*);

/*
__global__ void daxpy(int n, double a, double *x, double *y);

__global__ void ddscal(const int n, const double a, double *y);

__global__ void spmv_csr_half_vector_kernel(int, int*, int*, double*, double*, double*, int);

*/

void spmv_jad(matrix_t *mat, double *x, double *y, int neg);

void spmv_dia(matrix_t *mat, double *x, double *y, int neg);

void spmv_csr_cpu(csr_t *, double *, double*);

void dump_mat_coo(csr_t *A, char*);

void remove_diag(csr_t *);

void cuda_check_err();

void Partition1D(int len, int pnum, int idx, int &j1, int &j2);

void filtered_lanczos(matrix_t_*, int, int, double, double, int);

extern "C" void csrjad_(int *, double *, 
int *, int *, int *, int *, double *,
int *, int *);

extern "C" void multic_(int*, int*, int*, 
int*, int*, int*,int*, int*, int*);

extern "C" void aplb_(int*, int*, int*, 
double*, int*, int*,
double*, int*, int*, double*, int*, int*,
int*, int*, int*);

extern "C" void filter_(int*, int*, double*, 
double*, int*, int*, double*, 
int*, int*, int*, int*);

extern "C" void dperm_(int*, double *, int *, 
int*, double *,int *, int *, int *, int *, int*);

extern "C" void vperm_(int*, double *, int*);

extern "C" void getu_(int*,double*,int*,
int*,double*,int*,int*);

extern "C" void submat_(int*, int*, int*, 
int*, int*, int*, double*, int*, 
int*, int*, int*, double*, int*, int*);

extern "C" void csrdia_(int*, int*, int*, 
double*, int*, int*, int*, double *, 
int*, double*, int *, int*, int*);

extern "C" void METIS_PartGraphKway
(int*, int*, int*, int*, int*, 
int*, int*, int*, int*, int*, int*);

extern "C" void METIS_NodeND
(int*, int*, int*, int*, int*, int*, int*);

extern "C" void dsteqr_
(char*, int*, double*, double*, double*, int*, 
double*, int*);

extern "C" void daxpy_
(int *, double *, double *, int *, double *, int*);

/*----------------------------------------*/

#define STEQR dsteqr_
#define DAXPY daxpy_
#define GGEVX dggevx_
#define CUDOT cublasDdot
#define CUAXPY cublasDaxpy
#define CUSCAL cublasDscal
#define CUNRM2 cublasDnrm2
#define CUGEMV cublasDgemv
#define CUGEMM cublasDgemm

#define CHECKERR(ierr) assert(!(ierr))
#define Malloc(base, nmem, type) {\
  (base) = (type *)malloc((nmem)*sizeof(type)); \
  CHECKERR((base) == NULL); \
}
#define Calloc(base, nmem, type) {\
  (base) = (type *)calloc((nmem), sizeof(type)); \
  CHECKERR((base) == NULL); \
}
#define Realloc(base, nmem, type) {\
  (base) = (type *)realloc((base), (nmem)*sizeof(type)); \
  CHECKERR((base) == NULL && nmem > 0); \
}

#define CUDA_SAFE_CALL_NO_SYNC( call) {                               \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
    } }
#define CUDA_SAFE_CALL( call)    CUDA_SAFE_CALL_NO_SYNC(call);

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

