include ../../make.inc

COBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu))

all: $(COBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(INCS)

clean:
	-rm *.o
	
	
	
	
