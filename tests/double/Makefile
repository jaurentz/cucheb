HOMEDIR := /home/jared/Dropbox/MySoftware
CUCHEBDIR := $(HOMEDIR)/cucheb/cucheb.0.1
CUDADIR := /usr/local/cuda
CULADIR := /usr/local/cula
CULASPARSEDIR := /usr/local/culasparse

CFLAGS := -fopenmp

CUC := $(CUDADIR)/bin/nvcc
CUFLAGS := -O3 --compiler-options '$(CFLAGS)'\
		   -gencode arch=compute_20,code=sm_20\
		   -gencode arch=compute_30,code=sm_30\
		   -gencode arch=compute_35,code=sm_35\
		   #-gencode arch=compute_13,code=sm_13

INCLUDES := -I $(CUCHEBDIR)/include -I $(CUDADIR)/include -I $(CULADIR)/include #-I $(CULASPARSEDIR)/include
DIRS :=  -L $(CUCHEBDIR)/lib -L $(CULADIR)/lib64
LIBS := -lcucheb -lcuda -lcudart -lcublas -lcufft -lcurand -lcusparse -llapack -lblas -lm -lcula_lapack

BINS := $(patsubst %.cu,%,$(wildcard *.cu))

all: $(BINS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(INCLUDES) $(DIRS) $(LIBS)
	
clean:
	-rm $(BINS)
