include ../../environment

all: $(COBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< $(INCDIRS) $(LIBDIRS) $(LIBS)

clean:
	-rm *.o
	
	
	
	
