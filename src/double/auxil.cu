#include <gpusollib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#define MAXLINE 500

double wall_timer() {
   struct timeval tim;
   gettimeofday(&tim, NULL);
   double t = tim.tv_sec + tim.tv_usec/1e6;
   return(t);
}


/*----------------------------------------*/
void read_input( char* fn ) {

   FILE   *fp;
   char   s1[50],s2[50],s3[50],s4[50],s5[50];

/*----------------------------------------*/
  memset(s1, 0, 50*sizeof(char));
  memset(s2, 0, 50*sizeof(char));
  memset(s3, 0, 50*sizeof(char));
  memset(s4, 0, 50*sizeof(char));
  memset(s5, 0, 50*sizeof(char));

/*----------------------------------------*/
  fp = fopen(fn, "r");
  if (fp == NULL) {
     printf("Error in opening file %s\n",fn);
     exit(-1);
  }

  /*
  // matrix type
  fgets(line, MAXLINE, fp);
  assert(line[0] == '#');
  fgets(line, MAXLINE, fp);
  sscanf(line, "%s", s1);
  if (strcmp(s1, "csr") == 0)
    opts->mattype = CSR;
  else if (strcmp(s1, "jad") == 0)
    opts->mattype = JAD;
  else
    opts->mattype = DIA;
  */

  fclose(fp);
}


// Just print matrix
void dump_mat_coo(csr_t *A, char *fn) {
  int i,j,n;
  FILE *fp;
  fp = fopen(fn, "w");
  n = A->n;
  for (i=0; i<n; i++)
    for (j=(A->ia)[i]; j<(A->ia)[i+1]; j++)
      fprintf(fp, "%d %d %f\n", i+1, (A->ja)[j-1], (A->a)[j-1]);
  fclose(fp);
}










