# directory where CUCHEB will be installed
INSTALLDIR := $(HOME)/CUCHEB_INSTALL

# CUDA compiler and flags
CUC := /LUSTRE/apps/cuda-7.5/bin/nvcc
CUFLAGS := -O3 -gencode arch=compute_35,code=sm_35 

# CUDA libraries
CUDA := -lcudart 
CUBLAS := -lcublas 
CUSPARSE := -lcusparse
CUFFT := -lcufft 

# linking directories
INCS := 

# All libraries
LIBS := $(CUDA) $(CUBLAS) $(CUSPARSE) $(CUFFT)
