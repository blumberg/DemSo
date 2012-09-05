CC		:= nvcc
EXECUTABLE	:= DemSo
CUFILE		:= src/main.cu
CUHEADERS	:= src/main.cuh #main.h
CUFLAGS		:= -lglut -lpthread -arch=sm_11 -lGL
CUDEPS		:= src/functions.cu src/particles_kernel.cu

$(EXECUTABLE): $(CUFILE) $(CUDEPS) $(CUHEADERS)
	$(CC) $(CUFILE) -o $(EXECUTABLE) $(CUFLAGS) -include $(CUHEADERS)
	
clean:
	ls src/*~ | xargs rm -f
	ls doc/*~ | xargs rm -f
	ls includes/*~ | xargs rm -f
	ls *~ | xargs rm -f
	rm -f $(EXECUTABLE)
	
clean_tio:
	ls src/*~ | xargs rm -f
	ls includes/*~ | xargs rm -f
	ls *~ | xargs rm -f
