include ../make.inc

UTILDIR := ./utilities
SSRCDIR := ./single
CSRCDIR := ./complex_single
DSRCDIR := ./double
ZSRCDIR := ./complex_double

UOBJS := $(patsubst $(UTILDIR)/%.cu,$(UTILDIR)/%.o,$(wildcard $(UTILDIR)/*.cu))
SOBJS := $(patsubst $(SSRCDIR)/%.cu,$(SSRCDIR)/%.o,$(wildcard $(SSRCDIR)/*.cu))
COBJS := $(patsubst $(CSRCDIR)/%.cu,$(CSRCDIR)/%.o,$(wildcard $(CSRCDIR)/*.cu))
DOBJS := $(patsubst $(DSRCDIR)/%.cu,$(DSRCDIR)/%.o,$(wildcard $(DSRCDIR)/*.cu))
ZOBJS := $(patsubst $(ZSRCDIR)/%.cu,$(ZSRCDIR)/%.o,$(wildcard $(ZSRCDIR)/*.cu))
OBJS := $(UOBJS) $(SOBJS) $(COBJS) $(DOBJS) $(ZOBJS)

# CUCHEB individual sources
USRCS := $(wildcard $(UTILDIR)/*.cu)
SSRCS := $(wildcard $(SSRCDIR)/*.cu)
CSRCS := $(wildcard $(CSRCDIR)/*.cu)
DSRCS := $(wildcard $(DSRCDIR)/*.cu)
ZSRCS := $(wildcard $(ZSRCDIR)/*.cu) 

all: $(OBJS)
	
$(UOBJS): $(USRCS)
	make -C $(UTILDIR)
	
$(USRCS):

$(SOBJS): $(SSRCS)
	make -C $(SSRCDIR)
	
$(SSRCS):

$(COBJS): $(CSRCS)
	make -C $(CSRCDIR)
	
$(CSRCS):

$(DOBJS): $(DSRCS)
	make -C $(DSRCDIR)
	
$(DSRCS):

$(ZOBJS): $(ZSRCS)
	make -C $(ZSRCDIR)
	
$(ZSRCS):

clean:
	make clean -C $(UTILDIR) &&\
	make clean -C $(SSRCDIR) &&\
	make clean -C $(CSRCDIR) &&\
	make clean -C $(DSRCDIR) &&\
	make clean -C $(ZSRCDIR)

