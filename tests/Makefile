include ../make.inc

STESTDIR := ./single
CTESTDIR := ./complex_single
DTESTDIR := ./double
ZTESTDIR := ./complex_double

STESTS := $(patsubst $(STESTDIR)/%.cu,$(STESTDIR)/%,$(wildcard $(STESTDIR)/*.cu))
CTESTS := $(patsubst $(CTESTDIR)/%.cu,$(CTESTDIR)/%,$(wildcard $(CTESTDIR)/*.cu))
DTESTS := $(patsubst $(DTESTDIR)/%.cu,$(DTESTDIR)/%,$(wildcard $(DTESTDIR)/*.cu))
ZTESTS := $(patsubst $(ZTESTDIR)/%.cu,$(ZTESTDIR)/%,$(wildcard $(ZTESTDIR)/*.cu))

# CUCHEB individual test sources
STESTSRCS := $(wildcard $(STESTDIR)/*.cu)
CTESTSRCS := $(wildcard $(CTESTDIR)/*.cu)
DTESTSRCS := $(wildcard $(DTESTDIR)/*.cu)
ZTESTSRCS := $(wildcard $(ZTESTDIR)/*.cu)  

TESTS := $(DTESTS)

all: $(TESTS)

$(STESTS): $(STESTSRCS)
	make -C $(STESTDIR)
	
$(STESTSRCS):

$(CTESTS): $(CTESTSRCS)
	make -C $(CTESTDIR)
	
$(CTESTSRCS):

$(DTESTS): $(DTESTSRCS)
	make -C $(DTESTDIR)
	
$(DTESTSRCS):

$(ZTESTS): $(ZTESTSRCS)
	make -C $(ZTESTDIR)

$(ZTESTSRCS):

clean:
	make clean -C $(STESTDIR) &&\
	make clean -C $(DTESTDIR) &&\
	make clean -C $(CTESTDIR) &&\
	make clean -C $(ZTESTDIR)
	
	
	
	
