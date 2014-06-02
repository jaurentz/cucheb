include ../make.inc

#all: $(STESTS) $(DTESTS) $(CTESTS) $(ZTESTS)
all: $(DTESTS)

#$(STESTS): $(STESTSRCS)
#	make -C $(STESTDIR)
	
#$(STESTSRCS):
 
$(DTESTS): $(DTESTSRCS)
	make -C $(DTESTDIR)
	
$(DTESTSRCS):

#$(CTESTS): $(CTESTSRCS)
#	make -C $(CTESTDIR)
	
#$(CTESTSRCS):

#$(ZTESTS): $(ZTESTSRCS)
#	make -C $(ZTESTDIR)
	
#$(ZTESTSRCS):

clean:
	make clean -C $(STESTDIR) &&\
	make clean -C $(DTESTDIR) &&\
	make clean -C $(CTESTDIR) &&\
	make clean -C $(ZTESTDIR)
	
	
	
	
