include ./environment

LIBNAME := cucheb
MAJOR := 0
MINOR := 1
VERSION := $(MAJOR).$(MINOR)

OBJS := $(UOBJS) $(SOBJS) $(DOBJS) $(COBJS) $(ZOBJS)
SRCS := $(USRCS) $(SSRCS) $(DSRCS) $(CSRCS) $(ZSRCS)
HDRS := $(wildcard $(CUCHEBDIR)/include/*.h)

all: lib$(LIBNAME).so.$(VERSION)

install: lib$(LIBNAME).so.$(VERSION)
	mkdir -p $(INSTALLDIR)/cucheb 
	mkdir -p $(INSTALLDIR)/cucheb/include
	mkdir -p $(INSTALLDIR)/cucheb/lib
	cp -r $(CUCHEBDIR)/include/*.h $(INSTALLDIR)/cucheb/include/
	mv $(CUCHEBDIR)/lib$(LIBNAME).so.$(VERSION) $(INSTALLDIR)/cucheb/lib/
	ln -s $(INSTALLDIR)/cucheb/lib/lib$(LIBNAME).so.$(VERSION) lib$(LIBNAME).so
	mv $(CUCHEBDIR)/lib$(LIBNAME).so $(INSTALLDIR)/cucheb/lib/
	touch /etc/ld.so.conf.d/cucheb.conf
	echo "$(INSTALLDIR)/cucheb/lib" > /etc/ld.so.conf.d/cucheb.conf
	ldconfig

uninstall:
	-rm -rf $(INSTALLDIR)/cucheb /etc/ld.so.conf.d/cucheb.conf
 
lib$(LIBNAME).so.$(VERSION): $(OBJS)
	$(CUC) $(SOFLAGS) -o lib$(LIBNAME).so.$(VERSION) $^

$(OBJS): $(HDRS) $(SRCS)
	make -C $(CUCHEBDIR)/src

$(HDRS):

$(SRCS):

clean:
	-rm $(CUCHEBDIR)/lib$(LIBNAME).so.$(VERSION) &&\
	make clean -C $(CUCHEBDIR)/src
	
