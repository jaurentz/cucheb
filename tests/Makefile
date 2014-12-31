include ../make.inc

DTESTDIR := ./double

DTESTS := $(patsubst $(DTESTDIR)/%.cu,$(DTESTDIR)/%,$(wildcard $(DTESTDIR)/*.cu))

# CUCHEB individual test sources
DTESTSRCS := $(wildcard $(DTESTDIR)/*.cu)

TESTS := $(DTESTS)

all: $(TESTS)

$(DTESTS): $(DTESTSRCS)
	@$(MAKE) -C $(DTESTDIR)
	
$(DTESTSRCS):

clean:
	@$(MAKE) clean -C $(DTESTDIR)
	
	
	
	
