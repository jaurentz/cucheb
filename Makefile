include ./make.inc

# cucheb lib name
LIBNAME := cucheb
MAJOR := 0
MINOR := 1
PATCH := 0
VERSION := $(MAJOR).$(MINOR).$(PATCH)
export LIBNAME
export VERSION

# CUCHEB individual sources
SRCS := $(wildcard ./src/*/*.cu) $(wildcard ./src/*/*.f)

# CUCHEB individual objects
OBJS := $(SRCS:.cu=.o)
OBJS := $(OBJS:.f=.o)

# CUCHEB individual tests
TESTSRCS := $(wildcard ./tests/*/*.cu)
TESTS := $(TESTSRCS:.cu=)

all: lib$(LIBNAME).so.$(VERSION)

tests: $(TESTS) $(TESTSRCS)
	@$(MAKE) -C ./tests
	
$(TESTSRCS):

$(TESTS):

lib$(LIBNAME).so.$(VERSION): $(OBJS)
	@$(CUC) --compiler-options '-fPIC' --shared -o $@ $^

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
	@mv ./lib$(LIBNAME).so $(INSTALLDIR)/cucheb/lib/
	
uninstall: clean
	@rm -rf $(INSTALLDIR)/cucheb

clean:
	@$(MAKE) clean -C ./src
	@$(MAKE) clean -C ./tests

	
