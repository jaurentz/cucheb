include ../make.inc

DTESTSRCS := $(wildcard ./double/*.cu)
DTESTS := $(DTESTSRCS:.cu=)

all: $(DTESTS)

$(DTESTS): $(DTESTSRCS)
	@$(MAKE) -C ./double
	
$(DTESTSRCS):

clean:
	@$(MAKE) clean -C ./double
	
	
	
	
