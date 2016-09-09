include ../../make.inc

EXS := $(wildcard ./*.cu)
EXS := $(EXS:.cu=)

CUCHEBINCS := -I $(INSTALLDIR)/cucheb/include
CUCHEBLIBS := -L $(INSTALLDIR)/cucheb/lib -lcucheb 

all: $(EXS)

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(CUCHEBINCS) $(CUCHEBLIBS) $(INCS) $(LIBS) 
	
clean:
	@rm -f $(EXS)
