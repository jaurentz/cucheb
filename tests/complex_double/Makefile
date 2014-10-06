include ../../make.inc

ZTESTS := $(patsubst ./%.cu,./%,$(wildcard ./*.cu))

INSTALLEDINCLUDES := -I $(INSTALLDIR)/cucheb/include
INSTALLEDLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb

all: $(ZTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(INSTALLEDINCLUDES) $(INSTALLEDLIBS) $(CUDA) $(LAPACKE)
	
clean:
	rm $(ZTESTS)
