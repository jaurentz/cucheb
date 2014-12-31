# CUCHEB home directory
CUCHEBDIR := /usr1/jaurentz/Projects/cucheb

# directory where CUCHEB will be installed
INSTALLDIR := /usr1/jaurentz/MyLibs

# cucheb lib name
LIBNAME := cucheb
MAJOR := 0
MINOR := 1
VERSION := $(MAJOR).$(MINOR)

# C compiler and flags
CC := gcc
CFLAGS := -O3 -fPIC -fopenmp

# CUDA compiler and flags
CUC := nvcc
CUFLAGS := -O3 --compiler-options '$(CFLAGS)'\
			-gencode arch=compute_35,code=sm_35

# flags for building shared library
SOFLAGS := --compiler-options '-fPIC' --shared

# directory for CUDA library
CUDADIR := #/usr/local/cuda
CUDA := -lcuda -lcudart -lcublas -lcusparse -lcufft -lcurand 

# directory for LAPACKE library
LAPACKDIR := /usr1/jaurentz/MyLibs/lapack-3.5.0
LAPACK := 

# directory for LAPACKE library
LAPACKEDIR := $(LAPACKDIR)/lapacke
LAPACKE := -llapacke 

# linking directories
INCDIRS := -I $(CUCHEBDIR)/include -I $(LAPACKEDIR)/include 
INCS := $(INCDIRS) 
LIBDIRS :=   
LIBS := $(CUDA) -L $(LAPACKDIR) -L $(LAPACKEDIR) -llapacke -llapack -lblas -lm 

# CUCHEB src directory
SRCDIR := $(CUCHEBDIR)/src

# CUCHEB individual src directories
UDIR := $(SRCDIR)/utilities
DDIR := $(SRCDIR)/double

# CUCHEB individual srcs
USRCS := $(wildcard $(UDIR)/*.cu)
DSRCS := $(wildcard $(DDIR)/*.cu)

# CUCHEB individual objs
UOBJS := $(patsubst $(UDIR)/%.cu,$(UDIR)/%.o,$(USRCS))
DOBJS := $(patsubst $(DDIR)/%.cu,$(DDIR)/%.o,$(DSRCS))

# CUCHEB test directory
TESTDIR := $(CUCHEBDIR)/tests

# CUCHEB individual test directories
DTESTDIR := $(TESTDIR)/double

# CUCHEB individual test sources
DTESTSRCS := $(wildcard $(DTESTDIR)/*.cu)

# CUCHEB individual tests
DTESTS := $(patsubst $(DTESTDIR)/%.cu,%,$(DTESTSRCS))
