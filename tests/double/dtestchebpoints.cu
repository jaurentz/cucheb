#include <cucheb.h>

int main(){

	// compute variables
	/* double */
	int dn = pow(2,5)+1;
	int sizeD = sizeof(double);
	double *dpts;
	double da, db;
	
	// allocate host memory
	/* double */
	cuchebCheckError((void*)(dpts = (double*)malloc(dn*sizeD)),__FILE__,__LINE__);
	
	// set host pointers
	/* double */
	da = -1.0;
	db = 1.0;

	// compute chebpoints
	/* double */
	cuchebCheckError(cuchebDpoints(dn,&da,&db,dpts,1),__FILE__,__LINE__);
	
	// print output
	/* double */
	printf("\ndouble precision\n");
	for(int ii=0;ii<dn;ii++){printf("dpts[%d] = %+1.15e\n",ii,dpts[ii]);}
	printf("\n");

	// free host memory
	/* double */
	free(dpts);
	
	return 0;
}
