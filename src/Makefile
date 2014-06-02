include ../make.inc

#all: $(UOBJS) $(SOBJS) $(DOBJS) $(COBJS) $(ZOBJS)
all: $(UOBJS) $(DOBJS) 

$(UOBJS): $(USRCS)
	make -C $(UDIR)

$(USRCS):
 
#$(SOBJS): $(SSRCS)
#	make -C $(SDIR)

#$(SSRCS):

$(DOBJS): $(DSRCS)
	make -C $(DDIR)

$(DSRCS):

#$(COBJS): $(CSRCS)
#	make -C $(CDIR)

#$(CSRCS):

#$(ZOBJS): $(ZSRCS)
#	make -C $(ZDIR)

#$(ZSRCS):

clean:
	make clean -C $(UDIR) &&\
	make clean -C $(SDIR) &&\
	make clean -C $(DDIR) &&\
	make clean -C $(CDIR) &&\
	make clean -C $(ZDIR)
	
	
	
	
