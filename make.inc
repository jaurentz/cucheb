# directory where CUCHEB will be installed
INSTALLDIR := $(HOME)/Mylibs

# CUDA compiler and flags
CUC := /usr/local/cuda-7.0/bin/nvcc
CUFLAGS := -O3 -gencode arch=compute_35,code=sm_35 --std=c++11 

# CUDA libraries
CUDA := -lcuda 
CUBLAS := -lcublas 
CUSPARSE := -lcusparse
CUFFT := -lcufft 

# linking directories
INCS := 

# All libraries
LIBS := $(CUDA) $(CUBLAS) $(CUSPARSE) $(CUFFT)
