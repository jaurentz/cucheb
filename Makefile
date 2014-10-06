include ./make.inc

# cucheb lib name
LIBNAME := cucheb
MAJOR := 0
MINOR := 1
VERSION := $(MAJOR).$(MINOR)

# CUCHEB individual source directories
UTILDIR := ./src/utilities
SSRCDIR := ./src/single
CSRCDIR := ./src/complex_single
DSRCDIR := ./src/double
ZSRCDIR := ./src/complex_double

# CUCHEB individual objects
UOBJS := $(patsubst $(UTILDIR)/%.cu,$(UTILDIR)/%.o,$(wildcard $(UTILDIR)/*.cu))
SOBJS := $(patsubst $(SSRCDIR)/%.cu,$(SSRCDIR)/%.o,$(wildcard $(SSRCDIR)/*.cu))
COBJS := $(patsubst $(CSRCDIR)/%.cu,$(CSRCDIR)/%.o,$(wildcard $(CSRCDIR)/*.cu))
DOBJS := $(patsubst $(DSRCDIR)/%.cu,$(DSRCDIR)/%.o,$(wildcard $(DSRCDIR)/*.cu))
ZOBJS := $(patsubst $(ZSRCDIR)/%.cu,$(ZSRCDIR)/%.o,$(wildcard $(ZSRCDIR)/*.cu))
OBJS := $(UOBJS) $(DOBJS)

all: lib$(LIBNAME).so.$(VERSION)

tests: $(OBJS)
	make -C ./tests

lib$(LIBNAME).so.$(VERSION): $(OBJS)
	$(CUC) --compiler-options '-fPIC' --shared -o lib$(LIBNAME).so.$(VERSION) $^

$(OBJS):
	make -C ./src
	
install: lib$(LIBNAME).so.$(VERSION) 
	mkdir -p $(INSTALLDIR)/cucheb 
	mkdir -p $(INSTALLDIR)/cucheb/include
	mkdir -p $(INSTALLDIR)/cucheb/lib
	cp -r ./include/*.h $(INSTALLDIR)/cucheb/include/
	mv ./lib$(LIBNAME).so.$(VERSION) $(INSTALLDIR)/cucheb/lib/
	ln -s $(INSTALLDIR)/cucheb/lib/lib$(LIBNAME).so.$(VERSION) lib$(LIBNAME).so
	mv ./lib$(LIBNAME).so $(INSTALLDIR)/cucheb/lib/

uninstall: clean
	rm -rf $(INSTALLDIR)/cucheb

clean:
	make clean -C ./src &&\
	make clean -C ./tests

	
