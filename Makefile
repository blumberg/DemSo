CC		:= nvcc
EXECUTABLE	:= particles
CUFILE		:= src/main.cu
CUHEADERS	:= src/main.cuh #main.h
CUFLAGS		:= -lglut -lpthread
CUDEPS		:= src/functions.cu src/particles_kernel.cu

$(EXECUTABLE): $(CUFILE) $(CUDEPS) $(CUHEADERS)
	$(CC) $(CUFILE) -o $(EXECUTABLE) $(CUFLAGS) -include $(CUHEADERS)
	
clean:
	rm -f $(EXECUTABLE)
