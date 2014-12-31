include ../make.inc

UTILDIR := ./utilities
DSRCDIR := ./double

UOBJS := $(patsubst $(UTILDIR)/%.cu,$(UTILDIR)/%.o,$(wildcard $(UTILDIR)/*.cu))
DOBJS := $(patsubst $(DSRCDIR)/%.cu,$(DSRCDIR)/%.o,$(wildcard $(DSRCDIR)/*.cu))
OBJS := $(UOBJS) $(DOBJS)

# CUCHEB individual sources
USRCS := $(wildcard $(UTILDIR)/*.cu)
DSRCS := $(wildcard $(DSRCDIR)/*.cu)

all: $(OBJS)
	
$(UOBJS): $(USRCS)
	@$(MAKE) -C $(UTILDIR)
	
$(USRCS):

$(DOBJS): $(DSRCS)
	@$(MAKE) -C $(DSRCDIR)
	
$(DSRCS):

clean:
	@$(MAKE) clean -C $(UTILDIR)
	@$(MAKE) clean -C $(DSRCDIR)

