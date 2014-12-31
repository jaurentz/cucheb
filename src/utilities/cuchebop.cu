/*-----------------------------------------------------------------



-----------------------------------------------------------------*/

#include <cucheb.h>

/* constructors */
/* default */
ChebOp::ChebOp(void){
	
	// set field
	field = CUCHEB_FIELD_FLOAT;
	
	// set degree
	n = 1;
	
	// set opmult
	opmult = NULL;
	
	// set userdata
	userdata = NULL;
	
	// set chebpoly
	chebpoly = NULL;
}

/* user defined */
ChebOp::ChebOp(int N, cuchebOpMult OPMULT, void* USERDATA, ChebPoly* CHEBPOLY){

	// set n
	if(N > 0){n = N;}
	else{
		fprintf(stderr,"\nWarning in ChebOp constructor: n must be > 0!\n\n");
		cuchebExit(-1);
	}
	
	// set opmult
	if(OPMULT != NULL){opmult = OPMULT;}
	else{
		fprintf(stderr,"\nWarning in ChebOp constructor: opmult must != NULL!\n\n");
		cuchebExit(-1);
	}
	
	// set userdata
	userdata = USERDATA;
	
	// set chebpoly
	if(CHEBPOLY != NULL){chebpoly = CHEBPOLY;}
	else{
		fprintf(stderr,"\nWarning in ChebOp constructor: chebpoly must != NULL!\n\n");
		cuchebExit(-1);
	}
	
	// set field
	field = chebpoly->getField();
	if((field != CUCHEB_FIELD_FLOAT) && 
	(field != CUCHEB_FIELD_DOUBLE) &&
	(field != CUCHEB_FIELD_FLOAT_COMPLEX) &&
	(field != CUCHEB_FIELD_DOUBLE_COMPLEX)){
		fprintf(stderr,"\nWarning in ChebOp constructor: invalid field!\n\n");
		cuchebExit(-1);
	}	
}

/* copy */
ChebOp::ChebOp(const ChebOp& CHEBOP){

	// set n
	if(CHEBOP.getN() > 0){n = CHEBOP.getN();}
	else{
		fprintf(stderr,"\nWarning in ChebOp constructor: n must be > 0!\n\n");
		cuchebExit(-1);
	}
	
	// set opmult
	if(CHEBOP.getOpmult() != NULL){opmult = CHEBOP.getOpmult();}
	else{
		fprintf(stderr,"\nWarning in ChebOp constructor: opmult must != NULL!\n\n");
		cuchebExit(-1);
	}
	
	// set userdata
	userdata = CHEBOP.getUserdata();
	
	// set chebpoly
	if(CHEBOP.getChebpoly() != NULL){chebpoly = CHEBOP.getChebpoly();}
	else{
		fprintf(stderr,"\nWarning in ChebOp constructor: chebpoly must != NULL!\n\n");
		cuchebExit(-1);
	}
	
	// set field
	field = chebpoly->getField();
	if((field != CUCHEB_FIELD_FLOAT) && 
	(field != CUCHEB_FIELD_DOUBLE) &&
	(field != CUCHEB_FIELD_FLOAT_COMPLEX) &&
	(field != CUCHEB_FIELD_DOUBLE_COMPLEX)){
		fprintf(stderr,"\nWarning in ChebOp constructor: invalid field!\n\n");
		cuchebExit(-1);
	}	
}

/* mult */
cuchebStatus_t ChebOp::Mult(void*x, void* y){

	if(field == CUCHEB_FIELD_DOUBLE){
		cuchebCheckError(cuchebDmult(n,(double*)x,(double*)y,opmult,userdata,chebpoly),__FILE__,__LINE__);
	}
	else{
		fprintf(stderr,"\nWarning in ChebOp Mult: invalid field!\n\n");
		cuchebExit(-1);
	}
	
	return CUCHEB_STATUS_SUCCESS;
}

/* printers */
/* print short */
cuchebStatus_t ChebOp::print(void){

	// header
	printf("\nChebOp:\n");

	// field
	if(field == CUCHEB_FIELD_FLOAT){printf(" field = CUCHEB_FIELD_FLOAT\n");}
	else if(field == CUCHEB_FIELD_DOUBLE){printf(" field = CUCHEB_FIELD_DOUBLE\n");}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){printf(" field = CUCHEB_FIELD_FLOAT_COMPLEX\n");}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){printf(" field = CUCHEB_FIELD_DOUBLE_COMPLEX\n");}	

	// n
	printf(" n = %d",n);
	
	// chebpoly
	chebpoly->print();

	// footer
	printf("\n");
	
	return CUCHEB_STATUS_SUCCESS;
}

/* print long */
cuchebStatus_t ChebOp::printlong(void){

	// header
	printf("\nChebOp:\n");

	// field
	if(field == CUCHEB_FIELD_FLOAT){printf(" field = CUCHEB_FIELD_FLOAT\n");}
	else if(field == CUCHEB_FIELD_DOUBLE){printf(" field = CUCHEB_FIELD_DOUBLE\n");}
	else if(field == CUCHEB_FIELD_FLOAT_COMPLEX){printf(" field = CUCHEB_FIELD_FLOAT_COMPLEX\n");}
	else if(field == CUCHEB_FIELD_DOUBLE_COMPLEX){printf(" field = CUCHEB_FIELD_DOUBLE_COMPLEX\n");}	

	// n
	printf(" n = %d",n);
	
	// chebpoly
	chebpoly->printlong();

	// footer
	printf("\n");
	
	return CUCHEB_STATUS_SUCCESS;
}

// operator overload
ChebOp& ChebOp::operator= (const ChebOp& CHEBOP){

	// check for self-assignment
    if(this == &CHEBOP){return *this;}

	// set n
	if(CHEBOP.getN() > 0){n = CHEBOP.getN();}
	else{
		fprintf(stderr,"\nWarning in ChebOp assignment: n must be > 0!\n\n");
		cuchebExit(-1);
	}
	
	// set opmult
	if(CHEBOP.getOpmult() != NULL){opmult = CHEBOP.getOpmult();}
	else{
		fprintf(stderr,"\nWarning in ChebOp assignment: opmult must != NULL!\n\n");
		cuchebExit(-1);
	}
	
	// set userdata
	userdata = CHEBOP.getUserdata();
	
	// set chebpoly
	if(CHEBOP.getChebpoly() != NULL){chebpoly = CHEBOP.getChebpoly();}
	else{
		fprintf(stderr,"\nWarning in ChebOp assignment: chebpoly must != NULL!\n\n");
		cuchebExit(-1);
	}
	
	// set field
	field = chebpoly->getField();
	if((field != CUCHEB_FIELD_FLOAT) && 
	(field != CUCHEB_FIELD_DOUBLE) &&
	(field != CUCHEB_FIELD_FLOAT_COMPLEX) &&
	(field != CUCHEB_FIELD_DOUBLE_COMPLEX)){
		fprintf(stderr,"\nWarning in ChebOp assignment: invalid field!\n\n");
		cuchebExit(-1);
	}	
	
	// return
	return *this;
}

