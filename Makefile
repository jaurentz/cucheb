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
SRCS := $(wildcard ./src/*/*.cu)

# CUCHEB individual objects
OBJS := $(SRCS:.cu=.o)

all: lib$(LIBNAME).so.$(VERSION)

lib$(LIBNAME).so.$(VERSION): objects
	@$(CUC) --compiler-options '-fPIC' --shared -o $@ $(OBJS)

objects: FORCE
	@$(MAKE) -C ./src

tests: FORCE
	@$(MAKE) -C ./tests
	
FORCE:
	
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

	
