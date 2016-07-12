CC=mpif90
CFLAGS=-g -O0
milagro:minimilagro.f90
	$(CC) -o milagro $(CFLAGS) $<

run:milagro
	mpirun -np 4 milagro

clean:
	rm -f milagro *.o
