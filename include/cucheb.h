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
//#include <cula_lapack.h>
//#include <cula_lapack_device.h>
#include <lapacke.h>
using namespace std;
#include <cucheberror.h>
#include <cuchebpoly.h>
#include <cuchebop.h>
#include <cuchebsolve.h>
#include <cuchebeigs.h>

cuchebStatus_t cuchebSetGridBlocks(int n, dim3 *blockSize, dim3 *gridSize);

cuchebStatus_t cuchebDinit(int n,double *x,int incx,double val);

#endif /* __cucheb_h__ */
