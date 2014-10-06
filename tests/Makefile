include ../make.inc

STESTDIR := ./single
CTESTDIR := ./complex_single
DTESTDIR := ./double
ZTESTDIR := ./complex_double

STESTS := $(patsubst $(STESTDIR)/%.cu,$(STESTDIR)/%,$(wildcard $(STESTDIR)/*.cu))
CTESTS := $(patsubst $(CTESTDIR)/%.cu,$(CTESTDIR)/%,$(wildcard $(CTESTDIR)/*.cu))
DTESTS := $(patsubst $(DTESTDIR)/%.cu,$(DTESTDIR)/%,$(wildcard $(DTESTDIR)/*.cu))
ZTESTS := $(patsubst $(ZTESTDIR)/%.cu,$(ZTESTDIR)/%,$(wildcard $(ZTESTDIR)/*.cu))

TESTS := $(DTESTS)

all: $(TESTS)

$(STESTS):
	make -C $(STESTDIR)

$(CTESTS):
	make -C $(CTESTDIR)

$(DTESTS):
	make -C $(DTESTDIR)

$(ZTESTS):
	make -C $(ZTESTDIR)

clean:
	make clean -C $(STESTDIR) &&\
	make clean -C $(DTESTDIR) &&\
	make clean -C $(CTESTDIR) &&\
	make clean -C $(ZTESTDIR)
	
	
	
	
