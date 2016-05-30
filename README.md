# cucheb - CUDA accelerated large sparse eigensolvers #
Jared L. Aurentz, Vassilis Kalantzis and Yousef Saad, June 2016

## Introduction ##
__cucheb__ is a collection of C++ subroutines for accurately and 
efficiently solving large sparse matrix eigenvalue problems using 
NVIDIA brand GPUs. These methods are well suited for computing
eigenvalues of 2D/3D discretization matrices that arise in 
elliptic and parabolic PDEs.

### Current features ###
__cucheb-v0.1.0__ has the following features:
 - double precision eigensolvers for real symmetric matrices

## Installation ##
__cucheb__ is built on top of the [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
and a small number of C++ standard libraries. You must have the toolkit
installed before the library can be built.

### Linux ###
To install on a Linux machine simply move into the __cucheb__ root directory, 
edit the file __make.inc__ to suit your system and type:
```
make install
```
This creates a shared object library __libcucheb.so._version___ and copies 
it into the user specified installation directory. The installation does not 
create any symbolic links or export any library paths.

## Removing cucheb ##
If the source directory has not been removed simply move into the __cucheb__ 
root directory and type:
```
make uninstall
```
If the source directory has been removed the install directory will have to 
be removed explicitly by the user.

## Questions and issues ##
If you have any questions or encounter any issues while using __cucheb__ 
please file an issue on the [__cucheb__ issues](https://github.com/jaurentz
/cucheb/issues) page of Github.
