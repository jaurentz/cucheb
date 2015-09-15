include ../../make.inc

OBJS := $(wildcard ./*.cu)
OBJS := $(OBJS:.cu=.o)

all: $(OBJS)

%.o: %.cu $(INCS)
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC' -c $< -I ../../include $(INCS) 

$(INCS):

clean:
	@rm -f *.o
	
	
	
	
