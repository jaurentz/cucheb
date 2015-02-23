include ../../make.inc

UOBJS := $(wildcard ./*.cu)
UOBJS := $(UOBJS:.cu=.o)

all: $(UOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(INCS) 

clean:
	@rm -f *.o
	
	
	
	
