# directory where CUCHEB will be installed
INSTALLDIR := $(HOME)/Mylibs

# CUDA compiler and flags
CUC := /usr/local/cuda-7.0/bin/nvcc
CUFLAGS := -O3 -gencode arch=compute_35,code=sm_35 

# CUDA libraries
CUDA := -lcuda 
CUBLAS := -lcublas 
CUSPARSE := -lcusparse
CUFFT := -lcufft 
LAPACK := -L $(HOME)/Mylibs/lapack-3.5.0 -llapacke -llapack -lblas -lgfortran

# linking directories
INCS := -I/home/aurentz/Mylibs/lapack-3.5.0/lapacke/include 

# All libraries
LIBS := $(LAPACK) $(CUDA) $(CUBLAS) $(CUSPARSE) $(CUFFT)
