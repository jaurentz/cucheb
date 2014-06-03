/** \file cuchebpoly.h
 *  \brief Macros and prototypes for ChebPoly class.
 *
 *  This contains the macros and prototypes for the 
 *  ChebPoly class. 
 *
 *  \author Jared L. Aurentz
 *  \bug No known bugs.
 */
#ifndef __cuchebpoly_h__ 
#define __cuchebpoly_h__

#include <cucheberror.h>

/** \brief Maximum power of 2 for MAX_FLOAT_DEG.
 *  
 *  This is the largest power of 2 such that 
 *  MAX_FLOAT_DEG = 2^MAX_FLOAT_DEG_EXP.
 */
#ifdef MAX_FLOAT_DEG_EXP
#undef MAX_FLOAT_DEG_EXP
#define MAX_FLOAT_DEG_EXP 13
#else
#define MAX_FLOAT_DEG_EXP 13
#endif

/** \brief Maximum power of 2 for MAX_DOUBLE_DEG.
 *  
 *  This is the largest power of 2 such that 
 *  MAX_DOUBLE_DEG = 2^MAX_DOUBLE_DEG_EXP.
 */
#ifdef MAX_DOUBLE_DEG_EXP
#undef MAX_DOUBLE_DEG_EXP
#define MAX_DOUBLE_DEG_EXP 16
#else
#define MAX_DOUBLE_DEG_EXP 16
#endif

/** \brief Maximum degree for single precision.
 *  
 *  This is the maximum allowed degree for single 
 *  precision computations.
 */
#ifdef MAX_FLOAT_DEG
#undef MAX_FLOAT_DEG
#define MAX_FLOAT_DEG (int)pow(2,MAX_FLOAT_DEG_EXP)
#else
#define MAX_FLOAT_DEG (int)pow(2,MAX_FLOAT_DEG_EXP)
#endif

/** \brief Maximum degree for double precision.
 *  
 *  This is the maximum allowed degree for double 
 *  precision computations.
 */
#ifdef MAX_DOUBLE_DEG
#undef MAX_DOUBLE_DEG
#define MAX_DOUBLE_DEG (int)pow(2,MAX_DOUBLE_DEG_EXP)
#else
#define MAX_DOUBLE_DEG (int)pow(2,MAX_DOUBLE_DEG_EXP)
#endif

/** \brief Computational number fields.
 *  
 *  This is the number field in which the
 *  computations are performed.
 */
typedef enum cuchebField_t {
	/** \brief Real, single precision.
	 *  
	 *  This is equivalent to float in C/C++.
	 */
    CUCHEB_FIELD_FLOAT          = 0,
    /** \brief Real, double precision.
	 *  
	 *  This is equivalent to double in C/C++.
	 */
    CUCHEB_FIELD_DOUBLE         = 1,
	/** \brief Complex, single precision.
	 *  
	 *  This is equivalent to cuComplex in CUDA.
	 */
    CUCHEB_FIELD_FLOAT_COMPLEX  = 2,
	/** \brief Complex, double precision.
	 *  
	 *  This is equivalent to cuDoubleComplex in CUDA.
	 */
    CUCHEB_FIELD_DOUBLE_COMPLEX = 3
} cuchebField_t;

/** \brief Float functional form for interpolation.
 *  
 *  This is a function pointer to a function that computes \f$y[i*incy] = f(x[i*incx])\f$, \f$i=0,\ldots,n-1\f$,
 *  where x and y are pointers to float arrays allocated in device memory.
 *  @param n number of elements to be read from x and written to y, should only accept n > 0
 *  @param x pointer to float array allocated in device memory
 *  @param incx distance between elements in x, should only accept incx > 0
 *  @param y pointer to float array allocated in device memory
 *  @param incy distance between elements in y, should only accept incy > 0
 */
typedef cuchebStatus_t (*cuchebFloatFun) (int n,const float* x,int incx,float* y,int incy,void *userdata);

/** \brief Double functional form for interpolation.
 *  
 *  This is a function pointer to a function that computes \f$y[i*incy] = f(x[i*incx])\f$, \f$i=0,\ldots,n-1\f$,
 *  where x and y are pointers to double arrays allocated in device memory.
 *  @param n number of elements to be read from x and written to y, should only accept n > 0
 *  @param x pointer to double array allocated in device memory
 *  @param incx distance between elements in x, should only accept incx > 0
 *  @param y pointer to double array allocated in device memory
 *  @param incy distance between elements in y, should only accept incy > 0
 */
typedef cuchebStatus_t (*cuchebDoubleFun) (int n,const double* x,int incx,double* y,int incy,void *userdata);

/** \brief cuComplex functional form for interpolation.
 *  
 *  This is a function pointer to a function that computes \f$y[i*incy] = f(x[i*incx])\f$, \f$i=0,\ldots,n-1\f$,
 *  where x and y are pointers to cuComplex arrays allocated in device memory.
 *  @param n number of elements to be read from x and written to y, should only accept n > 0
 *  @param x pointer to cuComplex array allocated in device memory
 *  @param incx distance between elements in x, should only accept incx > 0
 *  @param y pointer to cuComplex array allocated in device memory
 *  @param incy distance between elements in y, should only accept incy > 0
 */
typedef cuchebStatus_t (*cuchebCuComplexFun) (int n,const cuComplex* x,int incx,cuComplex* y,int incy,void *userdata);

/** \brief cuDoubleComplex functional form for interpolation.
 *  
 *  This is a function pointer to a function that computes \f$y[i*incy] = f(x[i*incx])\f$, \f$i=0,\ldots,n-1\f$,
 *  where x and y are pointers to cuDoubleComplex arrays allocated in device memory.
 *  @param n number of elements to be read from x and written to y, should only accept n > 0
 *  @param x pointer to cuDoubleComplex array allocated in device memory
 *  @param incx distance between elements in x, should only accept incx > 0
 *  @param y pointer to cuDoubleComplex array allocated in device memory
 *  @param incy distance between elements in y, should only accept incy > 0
 */
typedef cuchebStatus_t (*cuchebCuDoubleComplexFun) (int n,const cuDoubleComplex* x,int incx,cuDoubleComplex* y,int incy,void *userdata);


/**
 *  \brief A class for Chebyshev polynomials. 
 *
 *  A class for constructing and storing Chebyshev polynomial interpolants
 *  represented in a basis of 
 *  <a href="http://en.wikipedia.org/wiki/Chebyshev_polynomials">Chebyshev polynomials of the first kind</a>.
 *  This class forms the backbone of CUCHEB.  
 */
class ChebPoly{

	private:
		/**
		 *  \brief Computational and storage field for ChebPoly. 
		 */
		cuchebField_t field;
		/**
		 *  \brief Degree of the ChebPoly.
		 */
		int degree;
		/**
		 *  \brief Domain of the ChebPoly.
		 *
		 *  This is a void pointer to device allocated memory.
		 *  The allocated memory will be of the type 
		 *  corresponding to \ref field. The values of a and b represent 
		 *  an interval \f$[a,b]\f$, \f$ a \neq b\f$.  
		 */
		void *a;
		/**
		 *  \brief Domain of the ChebPoly.
		 *
		 *  This is a void pointer to device allocated memory.
		 *  The allocated memory will be of the type 
		 *  corresponding to \ref field. The values of a and b represent 
		 *  an interval \f$[a,b]\f$, \f$ a \neq b\f$.  
		 */
		void *b;
		/**
		 *  \brief Coefficients of the ChebPoly.
		 *
		 *  This is a void pointer to device allocated memory.
		 *  The allocated memory will be of the type 
		 *  corresponding to \ref field and be of size \ref degree+1.
		 *  The coefficients are ordered from highest degree to lowest.
		 */
		void *coeffs;

	public:
		/**
		 *  \brief Default constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  corresponding to the constant function 1,
		 *  with \ref field = CUCHEB_FIELD_FLOAT.
		 */
		ChebPoly(void);
		/**
		 *  \brief Constant constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  corresponding to the constant function 1,
		 *  with \ref field = FIELD. 
		 *  @param FIELD must be a valid \ref cuchebField_t
		 */
		ChebPoly(cuchebField_t FIELD);
		/**
		 *  \brief Copy constructor for ChebPoly.
		 *
		 *  This constructor returns a copy of the ChebPoly 
		 *  referenced by CHEBPOLY.
		 *  @param CHEBPOLY can be any valid reference to a ChebPoly
		 */
		ChebPoly(const ChebPoly& CHEBPOLY);
		/**
		 *  \brief Fixed degree, float constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree = DEG, obtained by interpolating
		 *  FUN through the DEG+1 Chebyshev extreme points
		 *  in the interval [A,B], A != B. The \ref field
		 *  is set to CUCHEB_FIELD_FLOAT.
		 *  @param FUN user defined \ref cuchebFloatFun
		 *  @param A float pointer to host memory
		 *  @param B float pointer to host memory
		 *  @param DEG non-negative integer
		 */
		ChebPoly(cuchebFloatFun FUN,float* A,float* B,void *USERDATA,int DEG);
		/**
		 *  \brief User specified tolerance, float constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_FLOAT_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below the relative tolerance TOL > 0. The \ref field
		 *  is set to CUCHEB_FIELD_FLOAT.
		 *  @param FUN user defined \ref cuchebFloatFun
		 *  @param A float pointer to host memory
		 *  @param B float pointer to host memory
		 *  @param TOL float pointer to host memory
		 */
		ChebPoly(cuchebFloatFun FUN,float* A,float* B,void *USERDATA,float* TOL);
		/**
		 *  \brief Default tolerance, float constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_FLOAT_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below machine precision. The \ref field
		 *  is set to CUCHEB_FIELD_FLOAT.
		 *  @param FUN user defined \ref cuchebFloatFun
		 *  @param A float pointer to host memory
		 *  @param B float pointer to host memory
		 */
		ChebPoly(cuchebFloatFun FUN,float* A,float* B,void *USERDATA);
		/**
		 *  \brief Fixed degree, double constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree = DEG, obtained by interpolating
		 *  FUN through the DEG+1 Chebyshev extreme points
		 *  in the interval [A,B], A != B. The \ref field
		 *  is set to CUCHEB_FIELD_DOUBLE.
		 *  @param FUN user defined \ref cuchebDoubleFun
		 *  @param A double pointer to host memory
		 *  @param B double pointer to host memory
		 *  @param DEG non-negative integer
		 */
		ChebPoly(cuchebDoubleFun FUN,void *USERDATA,double* A,double* B,int DEG);
		/**
		 *  \brief User specified tolerance, double constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_DOUBLE_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below the relative tolerance TOL > 0. The \ref field
		 *  is set to CUCHEB_FIELD_DOUBLE.
		 *  @param FUN user defined \ref cuchebDoubleFun
		 *  @param A double pointer to host memory
		 *  @param B double pointer to host memory
		 *  @param TOL double pointer to host memory
		 */
		ChebPoly(cuchebDoubleFun FUN,void *USERDATA,double* A,double* B,double* TOL);
		/**
		 *  \brief Default tolerance, double constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_DOUBLE_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below machine precision. The \ref field
		 *  is set to CUCHEB_FIELD_DOUBLE.
		 *  @param FUN user defined \ref cuchebDoubleFun
		 *  @param A double pointer to host memory
		 *  @param B double pointer to host memory
		 */
		ChebPoly(cuchebDoubleFun FUN,void *USERDATA,double* A,double* B);
		/**
		 *  \brief Fixed degree, cuComplex constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree = DEG, obtained by interpolating
		 *  FUN through the DEG+1 Chebyshev extreme points
		 *  in the interval [A,B], A != B. The \ref field
		 *  is set to CUCHEB_FIELD_FLOAT_COMPLEX.
		 *  @param FUN user defined \ref cuchebCuComplexFun
		 *  @param A cuComplex pointer to host memory
		 *  @param B cuComplex pointer to host memory
		 *  @param DEG non-negative integer
		 */
		ChebPoly(cuchebCuComplexFun FUN,cuComplex* A,cuComplex* B,void *USERDATA,int DEG);
		/**
		 *  \brief User specified tolerance, cuComplex constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_FLOAT_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below the relative tolerance TOL > 0. The \ref field
		 *  is set to CUCHEB_FIELD_FLOAT_COMPLEX.
		 *  @param FUN user defined \ref cuchebCuComplexFun
		 *  @param A cuComplex pointer to host memory
		 *  @param B cuComplex pointer to host memory
		 *  @param TOL cuComplex pointer to host memory
		 */
		ChebPoly(cuchebCuComplexFun FUN,cuComplex* A,cuComplex* B,void *USERDATA,float* TOL);
		/**
		 *  \brief Default tolerance, cuComplex constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_FLOAT_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below machine precision. The \ref field
		 *  is set to CUCHEB_FIELD_FLOAT_COMPLEX.
		 *  @param FUN user defined \ref cuchebCuComplexFun
		 *  @param A cuComplex pointer to host memory
		 *  @param B cuComplex pointer to host memory
		 */
		ChebPoly(cuchebCuComplexFun FUN,cuComplex* A,cuComplex* B,void *USERDATA);
		/**
		 *  \brief Fixed degree, cuDoubleComplex constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree = DEG, obtained by interpolating
		 *  FUN through the DEG+1 Chebyshev extreme points
		 *  in the interval [A,B], A != B. The \ref field
		 *  is set to CUCHEB_FIELD_DOUBLE_COMPLEX.
		 *  @param FUN user defined \ref cuchebCuDoubleComplexFun
		 *  @param A cuDoubleComplex pointer to host memory
		 *  @param B cuDoubleComplex pointer to host memory
		 *  @param DEG non-negative integer
		 */
		ChebPoly(cuchebCuDoubleComplexFun FUN,cuDoubleComplex* A,cuDoubleComplex* B,void *USERDATA,int DEG);
		/**
		 *  \brief User specified tolerance, cuDoubleComplex constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_DOUBLE_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below the relative tolerance TOL > 0. The \ref field
		 *  is set to CUCHEB_FIELD_DOUBLE_COMPLEX.
		 *  @param FUN user defined \ref cuchebCuDoubleComplexFun
		 *  @param A cuDoubleComplex pointer to host memory
		 *  @param B cuDoubleComplex pointer to host memory
		 *  @param TOL cuDoubleComplex pointer to host memory
		 */
		ChebPoly(cuchebCuDoubleComplexFun FUN,cuDoubleComplex* A,cuDoubleComplex* B,void *USERDATA,double* TOL);
		/**
		 *  \brief Default tolerance, cuDoubleComplex constructor for ChebPoly.
		 *
		 *  This constructor returns the ChebPoly 
		 *  of \ref degree <= \ref MAX_DOUBLE_DEG, obtained by 
		 *  adaptively interpolating
		 *  FUN through a sequence Chebyshev grids
		 *  in the interval [A,B], A != B until the highest degree coefficients 
		 *  are below machine precision. The \ref field
		 *  is set to CUCHEB_FIELD_DOUBLE_COMPLEX.
		 *  @param FUN user defined \ref cuchebCuDoubleComplexFun
		 *  @param A cuDoubleComplex pointer to host memory
		 *  @param B cuDoubleComplex pointer to host memory
		 */
		ChebPoly(cuchebCuDoubleComplexFun FUN,cuDoubleComplex* A,cuDoubleComplex* B,void *USERDATA);
		/**
		 *  \brief Destructor for ChebPoly.
		 *
		 *  This destructor frees any memory allocated on the device.
		 */
		~ChebPoly(void);
		
		/**
		 *  \brief Basic printer for ChebPoly.
		 *
		 *  This prints out \ref field, \ref degree and [\ref a, \ref b].
		 */
		void print(void);
		/**
		 *  \brief Long printer for ChebPoly.
		 *
		 *  This prints out \ref field, \ref degree ,[\ref a, \ref b] 
		 *  and the entries of \ref coeffs.
		 */
		void printlong(void);
		
		/**
		 *  \brief Returns \ref field for a ChebPoly.
		 *
		 *  Returns \ref field for a ChebPoly.
		 */
		cuchebField_t getField(void) const{return field;}
		/**
		 *  \brief Returns \ref degree for a ChebPoly.
		 *
		 *  Returns \ref degree for a ChebPoly.
		 */
		int getDegree(void) const{return degree;}
		/**
		 *  \brief Returns \ref a for a ChebPoly.
		 *
		 *  Returns \ref a for a ChebPoly.
		 */
		void* getA(void) const{return a;}
		/**
		 *  \brief Returns \ref b for a ChebPoly.
		 *
		 *  Returns \ref b for a ChebPoly.
		 */
		void* getB(void) const{return b;}
		/**
		 *  \brief Returns \ref coeffs for a ChebPoly.
		 *
		 *  Returns \ref coeffs for a ChebPoly.
		 */
		void* getCoeffs(void) const{return coeffs;}
		
		/**
		 *  \brief Assignment operator overload for a ChebPoly.
		 *
		 *  Assignment operator overload for a ChebPoly.
		 */
		ChebPoly& operator= (const ChebPoly&);
};

/**
 *  \brief Single precision Chebyshev extreme points for the interval [a,b].
 *
 *  This computes the Chebyshev extreme points of type float
 *  for the interval [a,b].
 *  @param n number of extreme points, must be > 1
 *  @param a float pointer to device memory
 *  @param b float pointer to device memory
 *  @param pts float pointer to device memory where the extreme points will be stored
 *  @param incpts distance between elements in pts, must be > 0
 */
cuchebStatus_t cuchebSpoints (int n, 
							  const float *a, 
							  const float *b, 
							  float *pts, 
							  int incpts);
/**
 *  \brief Double precision Chebyshev extreme points for the interval [a,b].
 *
 *  This computes the Chebyshev extreme points of type double
 *  for the interval [a,b].
 *  @param n number of extreme points, must be > 1
 *  @param a double pointer to device memory
 *  @param b double pointer to device memory
 *  @param pts double pointer to device memory where the extreme points will be stored
 *  @param incpts distance between elements in pts, must be > 0
 */
cuchebStatus_t cuchebDpoints (int n,
							  const double *a, 
							  const double *b, 
							  double *pts, 
						      int incpts);
/**
 *  \brief Complex single precision Chebyshev extreme points for the interval [a,b].
 *
 *  This computes the Chebyshev extreme points of type cuComplex
 *  for the interval [a,b].
 *  @param n number of extreme points, must be > 1
 *  @param a cuComplex pointer to device memory
 *  @param b cuComplex pointer to device memory
 *  @param pts cuComplex pointer to device memory where the extreme points will be stored
 *  @param incpts distance between elements in pts, must be > 0
 */
cuchebStatus_t cuchebCpoints (int n, 
							  const cuComplex *a, 
						      const cuComplex *b, 
							  cuComplex *pts, 
			 			      int incpts);
/**
 *  \brief Complex double precision Chebyshev extreme points for the interval [a,b].
 *
 *  This computes the Chebyshev extreme points of type cuDoubleComplex
 *  for the interval [a,b].
 *  @param n number of extreme points, must be > 1
 *  @param a cuDoubleComplex pointer to device memory
 *  @param b cuDoubleComplex pointer to device memory
 *  @param pts cuDoubleComplex pointer to device memory where the extreme points will be stored
 *  @param incpts distance between elements in pts, must be > 0
 */
cuchebStatus_t cuchebZpoints (int n, 
							  const cuDoubleComplex *a, 
							  const cuDoubleComplex *b, 
							  cuDoubleComplex *pts, 
							  int incpts);
							  
/**
 *  \brief Single precision Chebyshev coefficients.
 *
 *  This computes the Chebyshev coefficients of type float for a polynomial 
 *  interpolant through the Chebyshev extreme points.
 *  @param n number of extreme points, must be > 1
 *  @param fvals float pointer to device memory where the function values at the extreme points are stored
 *  @param incfvals distance between elements in fvals, must be > 0
 *  @param coeffs float pointer to device memory where the coefficients will be stored
 *  @param inccfs distance between elements in coeffs, must be > 0
 */
cuchebStatus_t cuchebScoeffs (int n, 
							  const float *fvals,
							  int incfvals, 
							  float *coeffs, 
							  int inccfs);
/**
 *  \brief Double precision Chebyshev coefficients.
 *
 *  This computes the Chebyshev coefficients of type double for a polynomial 
 *  interpolant through the Chebyshev extreme points.
 *  @param n number of extreme points, must be > 1
 *  @param fvals double pointer to device memory where the function values at the extreme points are stored
 *  @param incfvals distance between elements in fvals, must be > 0
 *  @param coeffs double pointer to device memory where the coefficients will be stored
 *  @param inccfs distance between elements in coeffs, must be > 0
 */
cuchebStatus_t cuchebDcoeffs (int n,
							  const double *fvals, 
							  int incfvals, 
							  double *coeffs, 
						      int inccfs);
/**
 *  \brief Complex single precision Chebyshev coefficients.
 *
 *  This computes the Chebyshev coefficients of type cuComplex for a polynomial 
 *  interpolant through the Chebyshev extreme points.
 *  @param n number of extreme points, must be > 1
 *  @param fvals cuComplex pointer to device memory where the function values at the extreme points are stored
 *  @param incfvals distance between elements in fvals, must be > 0
 *  @param coeffs cuComplex pointer to device memory where the coefficients will be stored
 *  @param inccfs distance between elements in coeffs, must be > 0
 */
cuchebStatus_t cuchebCcoeffs (int n, 
							  const cuComplex *fvals, 
						      int incfvals, 
							  cuComplex *coeffs, 
			 			      int inccfs);
/**
 *  \brief Complex double precision Chebyshev coefficients.
 *
 *  This computes the Chebyshev coefficients of type cuDoubleComplex for a polynomial 
 *  interpolant through the Chebyshev extreme points.
 *  @param n number of extreme points, must be > 1
 *  @param fvals cuDoubleComplex pointer to device memory where the function values at the extreme points are stored
 *  @param incfvals distance between elements in fvals, must be > 0
 *  @param coeffs cuDoubleComplex pointer to device memory where the coefficients will be stored
 *  @param inccfs distance between elements in coeffs, must be > 0
 */
cuchebStatus_t cuchebZcoeffs (int n, 
							  const cuDoubleComplex *fvals, 
							  int incfvals, 
							  cuDoubleComplex *coeffs, 
							  int inccfs);


#endif /* __cuchebpoly_h__ */
