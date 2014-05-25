include ../../environment

all: $(SOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< $(INCDIRS) $(LIBDIRS) $(LIBS)

clean:
	-rm *.o
	
	
	
	
