include ../../make.inc

OBJS := $(wildcard ./*.cu)
OBJS := $(OBJS:.cu=.o)

all: $(OBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I ../../include $(INCS) 

clean:
	@rm -f *.o
	
	
	
	
