include ../make.inc

USRCS := $(wildcard ./utilities/*.cu)
DSRCS := $(wildcard ./double/*.cu)

UOBJS := $(USRCS:.cu=.o)
DOBJS := $(DSRCS:.cu=.o)

all: $(UOBJS) $(DOBJS)
	
$(UOBJS): $(USRCS)
	@$(MAKE) -C ./utilities
	
$(USRCS):

$(DOBJS): $(DSRCS)
	@$(MAKE) -C ./double
	
$(DSRCS):

clean:
	@$(MAKE) clean -C ./utilities
	@$(MAKE) clean -C ./double

