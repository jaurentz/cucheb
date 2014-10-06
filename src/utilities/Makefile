include ../../make.inc

UOBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu))

all: $(UOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -I../../include $(INCS) -c $< $(CUDA) $(LAPACKE)

clean:
	rm *.o
	
	
	
	
