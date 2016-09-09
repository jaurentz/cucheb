# directory where CUCHEB will be installed
INSTALLDIR := $(HOME)/CUCHEB_INSTALL

# CUDA compiler and flags
CUC := nvcc
CUFLAGS := -O3 -gencode arch=compute_35,code=sm_35 

# CUDA libraries
CUDA := -lcuda 
CUBLAS := -lcublas 
CUSPARSE := -lcusparse
CUFFT := -lcufft 

# linking directories
INCS := 

# All libraries
LIBS := $(CUDA) $(CUBLAS) $(CUSPARSE) $(CUFFT)
