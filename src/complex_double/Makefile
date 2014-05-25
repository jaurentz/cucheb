include ../../environment

all: $(ZOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< $(INCDIRS) $(LIBDIRS) $(LIBS)

clean:
	-rm *.o
	
	
	
	
