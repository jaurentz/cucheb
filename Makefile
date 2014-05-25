CUCHEB := /home/jared/Dropbox/MySoftware/cucheb
CUDADIR := /usr/local/cuda
CULADIR := /usr/local/cula
CULASPARSEDIR := /usr/local/culasparse

CC := gcc
CFLAGS := -O3 -fPIC -fopenmp

CUC := $(CUDADIR)/bin/nvcc
CUFLAGS := -O3 --compiler-options '-fPIC -fopenmp' -arch sm_20

SOFLAGS := --compiler-options '-fPIC' --shared

INCLUDES := -I $(CUCHEBDIR)/include -I $(CUCHEBDIR)/src -I $(CUDADIR)/include -I $(CULADIR)/include -I $(CULASPARSEDIR)/include
DIRS :=  -L $(CULADIR)/lib64 -L $(CULASPARSEDIR)/lib64
LIBS := -lcuda -lcublas -lcudart -llapack -lblas -lm -lcufft -lcurand -lcula_sparse -lcusparse

LIBNAME := cucheb
MAJOR := 0
MINOR := 1
VERSION := $(MAJOR).$(MINOR)

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(patsubst %.cu,%.o,$(wildcard *.cu))
HDRS := $(patsubst %.cpp,%.h,$(wildcard *.cpp)) $(patsubst %.cu,%.h,$(wildcard *.cu)) cucheb.h

all: lib$(LIBNAME).so

lib$(LIBNAME).so: lib$(LIBNAME).so.$(VERSION)
	ldconfig -v -n .
	ln -s lib$(LIBNAME).so.$(MAJOR) lib$(LIBNAME).so
 
lib$(LIBNAME).so.$(VERSION): $(OBJS)
	$(CUC) $(SOFLAGS) -o lib$(LIBNAME).so.$(MAJOR) $^

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< $(DIRS) $(INCLUDES) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< $(DIRS) $(INCLUDES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< $(DIRS) $(INCLUDES) $(LIBS)

$(OBJS): $(HDRS)

$(HDRS):

clean:
	-rm *.o *.so*
