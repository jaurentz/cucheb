include ../../make.inc

DOBJS := $(wildcard ./*.cu)
DOBJS := $(DOBJS:.cu=.o)

all: $(DOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(CUDAINCS) 

clean:
	@rm -f *.o
	
	
	
	
