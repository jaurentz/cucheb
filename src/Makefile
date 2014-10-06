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

OBJS := $(UOBJS) $(DOBJS)

all: $(OBJS)
	
$(UOBJS):
	make -C $(UTILDIR)

$(SOBJS):
	make -C $(SSRCDIR)

$(COBJS):
	make -C $(CSRCDIR)

$(DOBJS):
	make -C $(DSRCDIR)

$(ZOBJS):
	make -C $(ZSRCDIR)

clean:
	make clean -C $(UTILDIR) &&\
	make clean -C $(SSRCDIR) &&\
	make clean -C $(CSRCDIR) &&\
	make clean -C $(DSRCDIR) &&\
	make clean -C $(ZSRCDIR)

