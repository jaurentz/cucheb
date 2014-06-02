include ../../make.inc

INSTALLEDINCLUDES := -I $(INSTALLDIR)/cucheb/include
INSTALLEDLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb

all: $(ZTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(INSTALLEDINCLUDES) $(INSTALLEDLIBS) $(INCDIRS) $(LIBDIRS) $(LIBS)
	
clean:
	-rm $(ZTESTS)
