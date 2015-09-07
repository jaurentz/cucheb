include ../../make.inc

OBJS := $(wildcard ./*.cu)
OBJS := $(OBJS:.cu=.o)

all: $(OBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) --compiler-options '-fPIC' -c $< -I ../../include $(INCS) 

clean:
	@rm -f *.o
	
	
	
	
