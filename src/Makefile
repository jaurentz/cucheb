include ../environment

all: $(UOBJS) $(SOBJS) $(DOBJS) $(COBJS) $(ZOBJS)

$(UOBJS): $(USRCS)
	make -C $(UDIR)
 
$(SOBJS): $(SSRCS)
	make -C $(SDIR)

$(DOBJS): $(DSRCS)
	make -C $(DDIR)

$(COBJS): $(CSRCS)
	make -C $(CDIR)

$(ZOBJS): $(ZSRCS)
	make -C $(ZDIR)

clean:
	make clean -C $(UDIR) &&\
	make clean -C $(SDIR) &&\
	make clean -C $(DDIR) &&\
	make clean -C $(CDIR) &&\
	make clean -C $(ZDIR)
	
	
	
	
