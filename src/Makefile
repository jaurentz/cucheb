include ../make.inc

DSRCS := $(wildcard ./double/*.cu) $(wildcard ./double/*.f)
DOBJS := $(DSRCS:.cu=.o) 
DOBJS := $(DOBJS:.f=.o)

all: $(DOBJS)
	
$(DOBJS): $(DSRCS)
	@$(MAKE) -C ./double
	
$(DSRCS):

clean:
	@$(MAKE) clean -C ./double

