include ../../make.inc

ZOBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu))

all: $(ZOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(INCS)

clean:
	-rm *.o
	
	
	
	
