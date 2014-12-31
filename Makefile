include ./make.inc

# cucheb lib name
LIBNAME := cucheb
MAJOR := 0
MINOR := 1
VERSION := $(MAJOR).$(MINOR)

# CUCHEB individual source directories
UTILDIR := ./src/utilities
DSRCDIR := ./src/double

# CUCHEB individual sources
USRCS := $(wildcard $(UTILDIR)/*.cu)
DSRCS := $(wildcard $(DSRCDIR)/*.cu)
SRCS := $(USRCS) $(DSRCS) 

# CUCHEB individual objects
UOBJS := $(patsubst $(UTILDIR)/%.cu,$(UTILDIR)/%.o,$(wildcard $(UTILDIR)/*.cu))
DOBJS := $(patsubst $(DSRCDIR)/%.cu,$(DSRCDIR)/%.o,$(wildcard $(DSRCDIR)/*.cu))
OBJS := $(UOBJS) $(DOBJS)

# CUCHEB individual test directories
DTESTDIR := ./tests/double

# CUCHEB individual tests
DTESTSRCS := $(wildcard $(DTESTDIR)/*.cu)
TESTSRCS := $(DTESTSRCS)

all: lib$(LIBNAME).so.$(VERSION)

tests: $(OBJS) $(TESTSRCS)
	@$(MAKE) -C ./tests
	
$(TESTSRCS):

lib$(LIBNAME).so.$(VERSION): $(OBJS)
	$(CUC) --compiler-options '-fPIC' --shared -o lib$(LIBNAME).so.$(VERSION) $^

$(OBJS): $(SRCS)
	@$(MAKE) -C ./src

$(SRCS):
	
install: lib$(LIBNAME).so.$(VERSION) 
	@mkdir -p $(INSTALLDIR)/cucheb 
	@mkdir -p $(INSTALLDIR)/cucheb/include
	@mkdir -p $(INSTALLDIR)/cucheb/lib
	@cp -r ./include/*.h $(INSTALLDIR)/cucheb/include/
	@mv ./lib$(LIBNAME).so.$(VERSION) $(INSTALLDIR)/cucheb/lib/
	@ln -s $(INSTALLDIR)/cucheb/lib/lib$(LIBNAME).so.$(VERSION) lib$(LIBNAME).so
	mv ./lib$(LIBNAME).so $(INSTALLDIR)/cucheb/lib/

uninstall: clean
	@rm -rf $(INSTALLDIR)/cucheb

clean:
	@$(MAKE) clean -C ./src
	@$(MAKE) clean -C ./tests

	
