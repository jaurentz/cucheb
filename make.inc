# CUCHEB home directory
CUCHEBDIR := /home/jared/Projects/cuCheb

# directory where CUCHEB will be installed
INSTALLDIR := /usr/local

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
#CUC := /usr/local/cuda/bin/nvcc
CUFLAGS := -O3 --compiler-options '$(CFLAGS)'\
		   -gencode arch=compute_35,code=sm_35\
		   -gencode arch=compute_30,code=sm_30\
		   #-gencode arch=compute_20,code=sm_20
		   #-gencode arch=compute_13,code=sm_13

# flags for building shared library
SOFLAGS := --compiler-options '-fPIC' --shared

# directory for CUDA library
CUDADIR := #/usr/local/cuda
CUDA := -lcuda -lcudart -lcublas -lcusparse -lcufft -lcurand 

# directory for LAPACKE library
LAPACKDIR := 
LAPACK := #-llapack

# directory for LAPACKE library
LAPACKEDIR := /usr/local/lapacke
LAPACKE := -llapacke #-I $(LAPACKEDIR)/include -L $(LAPACKEDIR)/lib -llapacke

# linking directories
INCDIRS := -I $(CUCHEBDIR)/include #-I $(CUDADIR)/include -I $(LAPACKDIR)/include -I $(LAPACKEDIR)/include  
LIBDIRS :=   
LIBS := $(CUDA) $(LAPACK) $(LAPACKE) -lm 

# CUCHEB src directory
SRCDIR := $(CUCHEBDIR)/src

# CUCHEB individual src directories
UDIR := $(SRCDIR)/utilities
SDIR := $(SRCDIR)/single
DDIR := $(SRCDIR)/double
CDIR := $(SRCDIR)/complex_single
ZDIR := $(SRCDIR)/complex_double

# CUCHEB individual srcs
USRCS := $(wildcard $(UDIR)/*.cu)
SSRCS := $(wildcard $(SDIR)/*.cu)
DSRCS := $(wildcard $(DDIR)/*.cu)
CSRCS := $(wildcard $(CDIR)/*.cu)
ZSRCS := $(wildcard $(ZDIR)/*.cu)

# CUCHEB individual objs
UOBJS := $(patsubst $(UDIR)/%.cu,$(UDIR)/%.o,$(USRCS))
SOBJS := $(patsubst $(SDIR)/%.cu,$(SDIR)/%.o,$(SSRCS))
DOBJS := $(patsubst $(DDIR)/%.cu,$(DDIR)/%.o,$(DSRCS))
COBJS := $(patsubst $(CDIR)/%.cu,$(CDIR)/%.o,$(CSRCS))
ZOBJS := $(patsubst $(ZDIR)/%.cu,$(ZDIR)/%.o,$(ZSRCS))

# CUCHEB test directory
TESTDIR := $(CUCHEBDIR)/tests

# CUCHEB individual test directories
STESTDIR := $(TESTDIR)/single
DTESTDIR := $(TESTDIR)/double
CTESTDIR := $(TESTDIR)/complex_single
ZTESTDIR := $(TESTDIR)/complex_double

# CUCHEB individual test sources
STESTSRCS := $(wildcard $(STESTDIR)/*.cu)
DTESTSRCS := $(wildcard $(DTESTDIR)/*.cu)
CTESTSRCS := $(wildcard $(CTESTDIR)/*.cu)
ZTESTSRCS := $(wildcard $(ZTESTDIR)/*.cu)

# CUCHEB individual tests
STESTS := $(patsubst $(STESTDIR)/%.cu,%,$(STESTSRCS))
DTESTS := $(patsubst $(DTESTDIR)/%.cu,%,$(DTESTSRCS))
CTESTS := $(patsubst $(CTESTDIR)/%.cu,%,$(CTESTSRCS))
ZTESTS := $(patsubst $(ZTESTDIR)/%.cu,%,$(ZTESTSRCS))
