include ../../make.inc

DSRCS := $(wildcard ./*.cu) $(wildcard ./*.f)
DOBJS := $(DSRCS:.cu=.o) 
DOBJS := $(DOBJS:.f=.o)

all: $(DOBJS) 

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(CUDAINCS) 

%.o: %.f
	$(FC) $(FFLAGS) -fPIC -c $<

clean:
	@rm -f *.o
	
	
	
	
