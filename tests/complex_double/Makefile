include ../../make.inc

ZTESTS := $(patsubst ./%.cu,./%,$(wildcard ./*.cu))

INSTALLEDINCLUDES := -I$(INSTALLDIR)/cucheb/include
INSTALLEDLIBS := -L$(INSTALLDIR)/cucheb/lib -lcucheb

all: $(ZTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fopenmp' -o $@ $< $(INSTALLEDINCLUDES) $(INSTALLEDLIBS) $(INCS) $(CUDA) $(LAPACKE)
	
clean:
	-rm $(ZTESTS)
