# directory where CUCHEB will be installed
INSTALLDIR := /usr1/jaurentz/MyLibs

# CUDA compiler and flags
CUC := nvcc
CUFLAGS := -O3 -gencode arch=compute_35,code=sm_35

# CUDA libraries
CUDA := -lcuda -lcudart -lcublas -lcusparse -lcufft -lcurand 

# LAPACKE library
LAPACKE := -llapacke -llapack -lblas -lm

# linking directories
INCS := -I ./include -I $(LAPACKEDIR)/include 


