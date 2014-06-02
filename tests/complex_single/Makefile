include ../../make.inc

INSTALLEDINCLUDES := -I $(INSTALLDIR)/cucheb/include
INSTALLEDLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb

all: $(CTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(INSTALLEDINCLUDES) $(INSTALLEDLIBS) $(INCDIRS) $(LIBDIRS) $(LIBS)
	
clean:
	-rm $(CTESTS)
