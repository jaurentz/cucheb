include ../make.inc

TESTS := $(wildcard ./*.cu)
TESTS := $(TESTS:.cu=)

CUCHEBINCS := -I $(INSTALLDIR)/cucheb/include
CUCHEBLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb 

all: $(TESTS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(CUCHEBINCS) $(CUCHEBLIBS) $(INCS) $(LIBS) 
	
clean:
	@rm -f $(TESTS)
