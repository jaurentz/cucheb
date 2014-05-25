/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>

/* constructors */
/* default */
ChebPoly::ChebPoly(void){

	// compute variables
	float temp;
	
	// set field
	field = CUCHEB_FIELD_FLOAT;
	
	// set degree
	degree = 0;
	
	// set a and b
	cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
	cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
	temp = -1.0f;
	cuchebCheckError(cudaMemcpy(a, &temp, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	temp = 1.0f;
	cuchebCheckError(cudaMemcpy(b, &temp, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	
	// set coeffs
	cuchebCheckError(cudaMalloc(&coeffs, sizeof(float)),__FILE__,__LINE__);
	temp = 1.0f;
	cuchebCheckError(cudaMemcpy(coeffs, &temp, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
}

/* identity */
ChebPoly::ChebPoly(cuchebField_t F){

	// compute variables
	float stemp;
	double dtemp;
	cuComplex ctemp;
	cuDoubleComplex ztemp;
	
	// set field
	if(F == CUCHEB_FIELD_FLOAT){field = CUCHEB_FIELD_FLOAT;}
	else if(F == CUCHEB_FIELD_DOUBLE){field = CUCHEB_FIELD_DOUBLE;}
	else if(F == CUCHEB_FIELD_FLOAT_COMPLEX){field = CUCHEB_FIELD_FLOAT_COMPLEX;}
	else if(F == CUCHEB_FIELD_DOUBLE_COMPLEX){field = CUCHEB_FIELD_DOUBLE_COMPLEX;}
	else{printf("\nWarning in Chebpoly: Not a valid input for field! Setting field to CUCHEB_FIELD_FLOAT.\n\n");}
	
	// set a and b
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
		stemp = -1.0f;
		cuchebCheckError(cudaMemcpy(a, &stemp, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
		stemp = 1.0f;
		cuchebCheckError(cudaMemcpy(b, &stemp, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMalloc(&a, sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(double)),__FILE__,__LINE__);
		dtemp = -1.0;
		cuchebCheckError(cudaMemcpy(a, &dtemp, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
		dtemp = 1.0;
		cuchebCheckError(cudaMemcpy(b, &dtemp, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMalloc(&a, sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(cuComplex)),__FILE__,__LINE__);
		ctemp = make_cuFloatComplex(-1.0f,0.0f);
		cuchebCheckError(cudaMemcpy(a, &ctemp, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
		ctemp = make_cuFloatComplex(1.0f,0.0f);
		cuchebCheckError(cudaMemcpy(b, &ctemp, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMalloc(&a, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		ztemp = make_cuDoubleComplex(-1.0,0.0);
		cuchebCheckError(cudaMemcpy(a, &ztemp, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
		ztemp = make_cuDoubleComplex(1.0,0.0);
		cuchebCheckError(cudaMemcpy(b, &ztemp, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	
	// set coeffs
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(float)),__FILE__,__LINE__);
		stemp = 1.0f;
		cuchebCheckError(cudaMemcpy(coeffs, &stemp, sizeof(float), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(double)),__FILE__,__LINE__);
		dtemp = 1.0;
		cuchebCheckError(cudaMemcpy(coeffs, &dtemp, sizeof(double), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(cuComplex)),__FILE__,__LINE__);
		ctemp = make_cuFloatComplex(1.0f,0.0f);
		cuchebCheckError(cudaMemcpy(coeffs, &ctemp, sizeof(cuComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMalloc(&coeffs, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		ztemp = make_cuDoubleComplex(1.0,0.0);
		cuchebCheckError(cudaMemcpy(coeffs, &ztemp, sizeof(cuDoubleComplex), cudaMemcpyHostToDevice),__FILE__,__LINE__);
	}	
}

/* copy */
ChebPoly::ChebPoly(const ChebPoly& CP){

	// compute variables
	void *temp;

	// set field
	field = CP.getField();
	
	// set degree
	degree = CP.getDegree();
	
	// set a and b
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMalloc(&a, sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(double)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMalloc(&a, sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(cuComplex)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMalloc(&a, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// set coeffs
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
}

/* destructor */
ChebPoly::~ChebPoly(void){

	// free memory
	cuchebCheckError(cudaFree(a),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(b),__FILE__,__LINE__);
	cuchebCheckError(cudaFree(coeffs),__FILE__,__LINE__);
}

/* printers */
/* print short */
void ChebPoly::print(void){

	// compute variables
	float sa, sb;
	double da, db;
	cuComplex ca, cb;
	cuDoubleComplex za, zb;

	// header
	printf("\nChebPoly:\n");

	// field
	if(field == CUCHEB_FIELD_FLOAT){printf(" field = CUCHEB_FIELD_FLOAT\n");}
	else if(field == CUCHEB_FIELD_DOUBLE){printf(" field = CUCHEB_FIELD_DOUBLE\n");}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){printf(" field = CUCHEB_FIELD_FLOAT_COMPLEX\n");}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){printf(" field = CUCHEB_FIELD_DOUBLE_COMPLEX\n");}	

	// degree
	printf(" degree = %d\n",degree);
	
	// a and b
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMemcpy(&sa, a, sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&sb, b, sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [%+e,%+e]\n",sa,sb);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMemcpy(&da, a, sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&db, b, sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [%+1.15e,%+1.15e]\n",da,db);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMemcpy(&ca, a, sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&cb, b, sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [(%+e,%+e),(%+e,%+e)]\n",cuCrealf(ca),cuCimagf(ca),cuCrealf(cb),cuCimagf(cb));
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMemcpy(&za, a, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&zb, b, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [(%+1.15e,%+1.15e),(%+1.15e,%+1.15e)]\n",cuCreal(za),cuCimag(za),cuCreal(zb),cuCimag(zb));
	}

	// footer
	printf("\n");
}

/* print long */
void ChebPoly::printlong(void){

	// compute variables
	float sa, sb, *scoeffs;
	double da, db, *dcoeffs;
	cuComplex ca, cb, *ccoeffs;
	cuDoubleComplex za, zb, *zcoeffs;

	// header
	printf("\nChebPoly:\n");

	// field
	if(field == CUCHEB_FIELD_FLOAT){printf(" field = CUCHEB_FIELD_FLOAT\n");}
	else if(field == CUCHEB_FIELD_DOUBLE){printf(" field = CUCHEB_FIELD_DOUBLE\n");}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){printf(" field = CUCHEB_FIELD_FLOAT_COMPLEX\n");}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){printf(" field = CUCHEB_FIELD_DOUBLE_COMPLEX\n");}	

	// degree
	printf(" degree = %d\n",degree);
	
	// a and b
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMemcpy(&sa, a, sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&sb, b, sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [%+e,%+e]\n",sa,sb);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMemcpy(&da, a, sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&db, b, sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [%+1.15e,%+1.15e]\n",da,db);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMemcpy(&ca, a, sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&cb, b, sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [(%+e,%+e),(%+e,%+e)]\n",cuCrealf(ca),cuCimagf(ca),cuCrealf(cb),cuCimagf(cb));
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMemcpy(&za, a, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(&zb, b, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		printf(" [a,b] = [(%+1.15e,%+1.15e),(%+1.15e,%+1.15e)]\n",cuCreal(za),cuCimag(za),cuCreal(zb),cuCimag(zb));
	}
	
	// coeffs
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError((void*)(scoeffs = (float*)malloc((degree+1)*sizeof(float))),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(scoeffs, coeffs, (degree+1)*sizeof(float), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		for(int ii=0;ii<(degree+1);ii++){printf(" coeffs[%d] = %+e\n",ii,scoeffs[ii]);}
		free(scoeffs);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError((void*)(dcoeffs = (double*)malloc((degree+1)*sizeof(double))),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(dcoeffs, coeffs, (degree+1)*sizeof(double), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		for(int ii=0;ii<(degree+1);ii++){printf(" coeffs[%d] = %+1.15e\n",ii,dcoeffs[ii]);}
		free(dcoeffs);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError((void*)(ccoeffs = (cuComplex*)malloc((degree+1)*sizeof(cuComplex))),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(ccoeffs, coeffs, (degree+1)*sizeof(cuComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		for(int ii=0;ii<(degree+1);ii++){printf(" coeffs[%d] = (%+e,%+e)\n",ii,cuCrealf(ccoeffs[ii]),cuCimagf(ccoeffs[ii]));}
		free(ccoeffs);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError((void*)(zcoeffs = (cuDoubleComplex*)malloc((degree+1)*sizeof(cuDoubleComplex))),__FILE__,__LINE__);
		cuchebCheckError(cudaMemcpy(zcoeffs, coeffs, (degree+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost),__FILE__,__LINE__);
		for(int ii=0;ii<(degree+1);ii++){printf(" coeffs[%d] = (%+1.15e,%+1.15e)\n",ii,cuCreal(zcoeffs[ii]),cuCimag(zcoeffs[ii]));}
		free(zcoeffs);
	}

	// footer
	printf("\n");
}

// operator overload
ChebPoly& ChebPoly::operator= (const ChebPoly& CP){

	// check for self-assignment
    if(this == &CP){return *this;}

	// compute variables
	void *temp;

	// set field
	field = CP.getField();
	
	// set degree
	degree = CP.getDegree();
	
	// set a and b
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMalloc(&a, sizeof(float)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(float)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMalloc(&a, sizeof(double)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(double)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMalloc(&a, sizeof(cuComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(cuComplex)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMalloc(&a, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		cuchebCheckError(cudaMalloc(&b, sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		temp = CP.getA();
		cuchebCheckError(cudaMemcpy(a, temp, sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
		temp = CP.getB();
		cuchebCheckError(cudaMemcpy(b, temp, sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// set coeffs
	if(field == CUCHEB_FIELD_FLOAT){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(float)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(float), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(double)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(double), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuComplex)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(cuComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){
		cuchebCheckError(cudaMalloc(&coeffs, (degree+1)*sizeof(cuDoubleComplex)),__FILE__,__LINE__);
		temp = CP.getCoeffs();
		cuchebCheckError(cudaMemcpy(coeffs, temp, (degree+1)*sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice),__FILE__,__LINE__);
	}
	
	// return
	return *this;
}

