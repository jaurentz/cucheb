include ../make.inc

DSRCS := $(wildcard ./double/*.cu)

DOBJS := $(DSRCS:.cu=.o)

all: $(DOBJS)
	
$(DOBJS): $(DSRCS)
	@$(MAKE) -C ./double
	
$(DSRCS):

clean:
	@$(MAKE) clean -C ./double

