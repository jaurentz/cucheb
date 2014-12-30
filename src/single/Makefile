include ../../make.inc

SOBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu))

all: $(SOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(INCS)

clean:
	@rm -f *.o
	
	
	
	
