# directory where CUCHEB will be installed
INSTALLDIR := $(HOME)/Mylibs

# CUDA compiler and flags
CUC := nvcc #/usr/local/cuda-6.0/bin/nvcc
CUFLAGS := -O3 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_30,code=sm_30

# CUDA libraries
CUDADIR := #/usr/local/cuda-6.0
CUDA := -lcuda -lcudart -lcublas -lcusparse -lcufft -lcurand 

# LAPACKE library
LAPACKEDIR := $(HOME)/Mylibs/lapack-3.5.0
LAPACKE := -L$(LAPACKEDIR) -llapacke -llapack -lblas -lm

# linking directories
INCS := -I$(CUDADIR)/include -I$(LAPACKEDIR)/lapacke/include 


