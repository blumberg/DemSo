CC			:= nvcc
BINNAME		:= DemSo
CFLAGS		:= -arch=sm_11 -I../includes -DUSE_TEX
LIBS		:= -lglut -lpthread -lGL
OBJS		:= main.o functions.o particles_kernel.o

SUBDIRS 	:= src

export CC CFLAGS LIBS OBJS BINNAME

all:
	@for dir in $(SUBDIRS); do \
		(cd $$dir && $(MAKE) $(WHAT_TO_MAKE)) || exit 1; \
	done

#$(EXECUTABLE): $(CUFILE) $(CUDEPS) $(CUHEADERS)
#	$(CC) $(CUFILE) -o $(EXECUTABLE) $(CUFLAGS) -include $(CUHEADERS)

#debug: $(CUFILE) $(CUDEPS) $(CUHEADERS)
#	$(CC) -g -G $(CUFILE) -o $(EXECUTABLE) $(CUFLAGS) -include $(CUHEADERS)

#ptx: $(CUFILE) $(CUDEPS) $(CUHEADERS)
#	$(CC) $(CUFILE) -o $(EXECUTABLE).ptx $(CUFLAGS) -ptx -include $(CUHEADERS)

parser: src/parser.cpp src/DEMSimulation.cpp
	$(CC) src/parser.cpp src/DEMParticles.cpp src/DEMSimulation.cpp -o parser -I./includes $(CFLAGS)
	
clean:
	@for dir in $(SUBDIRS); do \
		(cd $$dir && $(MAKE) $@) || exit 1; \
	done
	
clean_tio:
	ls src/*~ | xargs rm -f
	ls includes/*~ | xargs rm -f
	ls *~ | xargs rm -f
