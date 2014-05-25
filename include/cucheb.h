/** \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *  
 * etc...
 */

#ifndef __cucheb_h__ /* __cucheb_h__ */
#define __cucheb_h__

#include <iostream>
#include <string>
#include <sstream>
#include <omp.h>
#include <cula_lapack.h>
#include <cula_lapack_device.h>
using namespace std;
#include <cucheberror.h>
#include <cuchebpoly.h>
#include <cuchebop.h>
#include <cuchebsolve.h>
#include <cuchebeig.h>

cuchebStatus_t cuchebSetGridBlocks(int n, dim3 *blockSize, dim3 *gridSize);

cuchebStatus_t cuchebSinit(int n,float *x,int incx,float val);
cuchebStatus_t cuchebDinit(int n,double *x,int incx,double val);
cuchebStatus_t cuchebCinit(int n,cuComplex *x,int incx,cuComplex val);
cuchebStatus_t cuchebZinit(int n,cuDoubleComplex *x,int incx,cuDoubleComplex val);

#endif /* __cucheb_h__ */