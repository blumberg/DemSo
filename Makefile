CC		:= nvcc
EXECUTABLE	:= DemSo
CUFILE		:= src/main.cu
CUHEADERS	:= src/main.cuh #main.h
CUFLAGS		:= -lglut -lpthread -arch=sm_11 -lGL -I./includes
CUDEPS		:= src/functions.cu src/particles_kernel.cu

$(EXECUTABLE): $(CUFILE) $(CUDEPS) $(CUHEADERS)
	$(CC) $(CUFILE) -o $(EXECUTABLE) $(CUFLAGS) -include $(CUHEADERS)

debug: $(CUFILE) $(CUDEPS) $(CUHEADERS)
	$(CC) -g -G $(CUFILE) -o $(EXECUTABLE) $(CUFLAGS) -include $(CUHEADERS)

ptx: $(CUFILE) $(CUDEPS) $(CUHEADERS)
	$(CC) $(CUFILE) -o $(EXECUTABLE).ptx $(CUFLAGS) -ptx -include $(CUHEADERS)

parser: src/parser.cpp src/DEMSimulation.cpp
	$(CC) src/parser.cpp src/DEMSimulation.cpp -o parser $(CUFLAGS)
	
clean:
	ls src/*~ | xargs rm -f
	ls doc/*~ | xargs rm -f
	ls includes/*~ | xargs rm -f
	ls *~ | xargs rm -f
	rm -f $(EXECUTABLE)
	rm -f $(EXECUTABLE).ptx
	
clean_tio:
	ls src/*~ | xargs rm -f
	ls includes/*~ | xargs rm -f
	ls *~ | xargs rm -f
