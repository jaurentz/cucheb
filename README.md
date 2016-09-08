# Cucheb - CUDA accelerated large sparse eigensolvers #
Jared L. Aurentz, Vassilis Kalantzis and Yousef Saad, September 2016

## Github ##
This README file is written in mark down. For the best experience please view this
file, along with the rest of the Cucheb library on
[github](https://github.com/jaurentz/cucheb).

## Introduction ##
Cucheb is a collection of C++ subroutines for accurately and efficiently
solving large sparse matrix eigenvalue problems using NVIDIA brand GPUs. These
methods are well suited for computing eigenvalues and eigenvectors of matrices
arising in quantum physics, structural engineering and network analysis.

### Current features ###
__cucheb-v0.1.0__ has the following features:
 - double precision eigensolvers for real symmetric matrices

## User-level programs ##
There are many files in the Cucheb library but only a few will be necessary for
most users. The first set of files are the objects used to store computed
quantities, such as eigenvalues and eigenvectors. Below we give a brief
description of each file and a link for further information:
 - [cuchebmatrix](https:/github.com/jaurentz/cucheb/src/double/cuchebmatrix_lanczos.cu)

## Installation ##
Cucheb is built on top of the [NVIDIA CUDA
Toolkit](https://developer.nvidia.com/cuda-toolkit) and a small number of C++
standard libraries. You must have the toolkit installed before the library can
be built.

### Linux ###
To install on a Linux machine simply move into the Cucheb root directory,
edit the file __make.inc__ to suit your system and type:
```
make install
```
This creates a shared object library __libcucheb.so._version___ and copies it
into the user specified installation directory. The installation does not
create any symbolic links or export any library paths.

## Removing Cucheb ##
If the source directory has not been removed simply move into the Cucheb
root directory and type:
```
make uninstall
```
If the source directory has been removed the install directory will have to be
removed explicitly by the user.

## Questions and issues ##
If you have any questions or encounter any issues while using Cucheb please
file an issue on the [Cucheb issues](https://github.com/jaurentz
/cucheb/issues) page of Github.
