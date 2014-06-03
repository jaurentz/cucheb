#include <cucheb.h>

/* chebpoints */
/* double precision */
cuchebStatus_t cuchebDpoints(int n,const double *a,const double *b,double *pts,int incpts){

	// check n
	if(n <= 1){
		fprintf(stderr,"\nIn %s line: %d, n must be > 1.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	if(n > MAX_DOUBLE_DEG+1){
		fprintf(stderr,"\nIn %s line: %d, n must be <= %d.\n",__FILE__,__LINE__,MAX_DOUBLE_DEG+1);
		cuchebExit(-1);
	}
	
	// check incpts
	if(incpts <= 0){
		fprintf(stderr,"\nIn %s line: %d, incpts must be > 0.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}

	// check a and b
	if(&a[0] == &b[0]){
		fprintf(stderr,"\nIn %s line: %d, a must not equal b.\n",__FILE__,__LINE__);
		cuchebExit(-1);
	}
	
	// chebpoints
	double theta;
	for(int ii=0;ii < n;ii++){
		theta = (double)(M_PI_2)*(2.0*ii-n+1.0)/(n-1.0);
		pts[ii*incpts] = (*b + *a)/2.0 + sin(theta)*(*b - *a)/2.0;
	}
	
	// return success
	return CUCHEB_STATUS_SUCCESS;
}
