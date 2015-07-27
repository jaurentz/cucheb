include ../../make.inc

DTESTS := $(wildcard ./*.cu)
DTESTS := $(DTESTS:.cu=)

CUCHEBINCS := -I $(INSTALLDIR)/cucheb/include
CUCHEBLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb -lcufft

all: $(DTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fopenmp' -o $@ $< $(CUCHEBINCS) $(CUCHEBLIBS) $(INCS) $(CUDALIBS) 
	
clean:
	@rm -f $(DTESTS)
