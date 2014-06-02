include ../../make.inc

INSTALLEDINCLUDES := -I $(INSTALLDIR)/cucheb/include
INSTALLEDLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb

all: $(DTESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< -L /usr/local/cucheb/lib -lcucheb $(INCDIRS) $(LIBDIRS) $(LIBS)
	
clean:
	-rm $(DTESTS)
