include ../../make.inc

DTESTS := $(wildcard ./*.cu)
DTESTS := $(DTESTS:.cu=)

CUCHEBINCS := -I $(INSTALLDIR)/cucheb/include
CUCHEBLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb -lcufft -lcublas -lcusparse -llapack -lblas

all: $(DTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fopenmp' -o $@ $< $(CUCHEBINCS) $(CUCHEBLIBS) $(CUDAINCS) $(CUDALIBS) -lgfortran -lm 
	
clean:
	@rm -f $(DTESTS)
