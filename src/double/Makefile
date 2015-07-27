include ../../make.inc

DOBJS := $(wildcard ./*.cu) $(wildcard ./*.f)
DOBJS := $(DOBJS:.cu=.o) $(DOBJS:.f=.o)

all: $(DOBJS) 

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC -fopenmp' -c $< -I../../include $(CUDAINCS) 

%.o: %.f
	$(FC) $(FFLAGS) -fPIC -c $<

clean:
	@rm -f *.o
	
	
	
	
