include ../make.inc

UTILDIR := ./utilities
SSRCDIR := ./single
DSRCDIR := ./double

UOBJS := $(patsubst $(UTILDIR)/%.cu,$(UTILDIR)/%.o,$(wildcard $(UTILDIR)/*.cu))
SOBJS := $(patsubst $(SSRCDIR)/%.cu,$(SSRCDIR)/%.o,$(wildcard $(SSRCDIR)/*.cu))
DOBJS := $(patsubst $(DSRCDIR)/%.cu,$(DSRCDIR)/%.o,$(wildcard $(DSRCDIR)/*.cu))
OBJS := $(UOBJS) $(SOBJS) $(DOBJS)

# CUCHEB individual sources
USRCS := $(wildcard $(UTILDIR)/*.cu)
SSRCS := $(wildcard $(SSRCDIR)/*.cu)
DSRCS := $(wildcard $(DSRCDIR)/*.cu)

all: $(OBJS)
	
$(UOBJS): $(USRCS)
	@$(MAKE) -C $(UTILDIR)
	
$(USRCS):

$(SOBJS): $(SSRCS)
	@$(MAKE) -C $(SSRCDIR)
	
$(SSRCS):

$(DOBJS): $(DSRCS)
	@$(MAKE) -C $(DSRCDIR)
	
$(DSRCS):

clean:
	@$(MAKE) clean -C $(UTILDIR)
	@$(MAKE) clean -C $(SSRCDIR)
	@$(MAKE) clean -C $(DSRCDIR)

