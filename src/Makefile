
DemSo: $(OBJS)
	$(CC) $(CFLAGS) -o $(BINNAME) $(OBJS) $(LIBS)

clean:
	-rm -f *~ *.o core $(BINNAME)

main.o: main.cu main.cuh
	$(CC) $(CFLAGS) -c -o $@ main.cu $(LIBS)

functions.o: functions.cu functions.cuh
	$(CC) $(CFLAGS) -c -o $@ functions.cu $(LIBS)

particles_kernel.o: particles_kernel.cu
	$(CC) $(CFLAGS) -c -o $@ particles_kernel.cu $(LIBS)
