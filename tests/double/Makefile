include ../../make.inc

DTESTS := $(patsubst ./%.cu,./%,$(wildcard ./*.cu))

INSTALLEDINCLUDES := -I$(INSTALLDIR)/cucheb/include
INSTALLEDLIBS := -L$(INSTALLDIR)/cucheb/lib -lcucheb

all: $(DTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fopenmp' -o $@ $< $(INSTALLEDINCLUDES) $(INSTALLEDLIBS) $(INCS) $(CUDA) $(LAPACKE)
	
clean:
	@rm -f $(DTESTS)
