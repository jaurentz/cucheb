include ../../make.inc

all: $(DOBJS)

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< $(INCDIRS) $(LIBDIRS) $(LIBS)

clean:
	-rm *.o
	
	
	
	
