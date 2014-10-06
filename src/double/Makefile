include ../../make.inc

DOBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu))

all: $(DOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -I../../include $(INCS) -c $< $(CUDA) $(LAPACKE)

clean:
	rm *.o
	
	
	
	
