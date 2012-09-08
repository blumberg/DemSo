
DemSo: $(OBJS)
	$(CC) $(CFLAGS) -o $(BINNAME) $(OBJS) $(LIBS)

clean:
	-rm -f *~ *.o core $(BINNAME)

main.o:
	$(CC) $(CFLAGS) -c -o $@ main.cu $(LIBS)

functions.o:
	$(CC) $(CFLAGS) -c -o $@ functions.cu $(LIBS)

particles_kernel.o:
	$(CC) $(CFLAGS) -c -o $@ particles_kernel.cu $(LIBS)
