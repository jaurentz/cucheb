include ../../make.inc

DOBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu))

all: $(DOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(INCS)

clean:
	-rm *.o
	
	
	
	
