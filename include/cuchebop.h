/** \file cuchebop.h
 *  \brief Macros and prototypes for ChebOp class.
 *
 *  This contains the macros and prototypes for the 
 *  ChebOp class. 
 *
 *  \author Jared L. Aurentz
 *  \bug No known bugs.
 */
#ifndef __cuchebop_h__ 
#define __cuchebop_h__

#include <cucheberror.h>

/** \brief General function pointer for matrix-vector multiplication.
 *  
 *  This is a function pointer to a function that computes y = Ax,
 *  for some user defined, square matrix A. The pointers x and y 
 *  are pointers to arrays allocated in device memory. These arrays can be
 *  allocated in any one of the types described by \ref cuchebField_t.
 *  @param x pointer to array allocated in device memory
 *  @param y pointer to array allocated in device memory
 *  @param userdata pointer to user defined data allocated in device memory
 */
typedef void (*cuchebOpMult)(void* x,void* y,void* userdata);

/**
 *  \brief A class for Chebyshev polynomials of Matrices. 
 *
 *  A class for storing and manipulating Chebyshev polynomials
 *  of large, sparse Matrices.
 */
class ChebOp{

	private:
		/**
		 *  \brief Computational field for ChebOp, required to cast the void pointers to the correct type.  
		 */
		cuchebField_t field;
		/**
		 *  \brief Dimension of Chebop, will be positive.  
		 */
		int n;
		/**
		 *  \brief Matrix-vector multiply routine to compute y = Ax.
		 *
		 *  The funciton that is pointed to by \ref opmult is not owned by ChebOp
		 *  and is neither allocated or freed via the constructor
		 *  or destructor.  
		 */
		cuchebOpMult opmult;
		/**
		 *  \brief User data that describes the matrix A.
		 *
		 *  The memory that is pointed to by \ref userdata is not owned by ChebOp
		 *  and is neither allocated or freed via the constructor
		 *  or destructor.  
		 */
		void* userdata;
		/**
		 *  \brief ChebPoly for computing y = p(A)x.
		 *
		 *  The memory that is pointed to by \ref chebpoly is not owned by ChebOp
		 *  and is neither allocated or freed via the constructor
		 *  or destructor.   
		 */
		ChebPoly* chebpoly;		
	
	public:
		/**
		 *  \brief Default constructor for ChebOp.
		 *
		 *  This constructor returns the ChebOp 
		 *  with \ref field = CUCHEB_FIELD_FLOAT,
		 *  \ref n = 1, \ref opmult = NULL, \ref userdata = 
		 *  NULL and \ref chebpoly = NULL.
		 */
		ChebOp(void);
		/**
		 *  \brief User data constructor for ChebOp.
		 *
		 *  This constructor returns the ChebOp 
		 *  with \ref n = N, \ref opmult = OPMULT,
		 *  \ref userdata = USERDATA, \ref chebpoly = 
		 *  CHEBPOLY and \ref field = CHEBPOLY->getField().
		 *  @param N dimension of ChebOp, must be > 0
		 *  @param OPMULT a valid \ref cuchebOpMult for computing y = Ax
		 *  @param USERDATA pointer to user defined data that describes A
		 *  @param CHEBPOLY pointer to a valid \ref ChebPoly object
		 */
		ChebOp(int N, cuchebOpMult OPMULT, void* USERDATA, ChebPoly* CHEBPOLY);
		/**
		 *  \brief Copy constructor for ChebOp.
		 *
		 *  This constructor returns a copy of the ChebOp 
		 *  that is referenced by CHEBOP.
		 *  @param CHEBOP reference to a valid ChebOp.
		 */
		ChebOp(const ChebOp& CHEBOP);
		
		/**
		 *  \brief Matrix-vector multiplier for ChebOp.
		 *
		 *  This computes y = p(A)x where y = Ax is
		 *  computed using \ref opmult.
		 *  @param x void pointer to device allocated memory of 
		 *  size \ref n and type \ref field.
		 *  @param y void pointer to device allocated memory of 
		 *  size \ref n and type \ref field.
		 */
		cuchebStatus_t Mult(void* x, void* y);
		
		/**
		 *  \brief Basic printer for ChebPoly.
		 *
		 *  This prints out \ref field, \ref n and \ref chebpoly.
		 */
		cuchebStatus_t print(void);
		/**
		 *  \brief Long printer for ChebPoly.
		 *
		 *  This prints out \ref field, \ref n and \ref chebpoly including 
		 *  the coefficients.
		 */
		cuchebStatus_t printlong(void);
		
		/**
		 *  \brief Returns \ref field for a ChebOp.
		 *
		 *  Returns \ref field for a ChebOp.
		 */
		cuchebField_t getField(void) const{return field;}
		/**
		 *  \brief Returns \ref n for a ChebOp.
		 *
		 *  Returns \ref n for a ChebOp.
		 */
		int getN(void) const{return n;}
		/**
		 *  \brief Returns \ref opmult for a ChebOp.
		 *
		 *  Returns \ref opmult for a ChebOp.
		 */
		cuchebOpMult getOpmult(void) const{return opmult;}
		/**
		 *  \brief Returns \ref userdata for a ChebOp.
		 *
		 *  Returns \ref userdata for a ChebOp.
		 */
		void* getUserdata(void) const{return userdata;}
		/**
		 *  \brief Returns \ref chebpoly for a ChebOp.
		 *
		 *  Returns \ref chebpoly for a ChebOp.
		 */
		ChebPoly* getChebpoly(void) const{return chebpoly;}
		
		/**
		 *  \brief Assignment operator overload for a ChebOp.
		 *
		 *  Assignment operator overload for a ChebOp.
		 */
		ChebOp& operator= (const ChebOp&);
};

/**
 *  \brief Double precision matrix-vector multiplier for ChebOp.
 *
 *  This computes y = p(A)x where A, x and y are
 *  stored as doubles. The product y = Ax is computed
 *  using opmult and userdata. 
 *  @param n dimension of the square matrix A
 *  @param x double pointer to device allocated memory
 *  @param y double pointer to device allocated memory 
 *  @param opmult function pointer to matrix-vector multiply routine
 *  @param userdata void pointer to memory that describes A
 *  @param chebpoly ChebPoly pointer to valid ChebPoly object
 */
cuchebStatus_t cuchebDmult(int n,double* x,double* y,cuchebOpMult opmult,void* userdata,ChebPoly* chebpoly);

#endif /* __cuchebop_h__ */
