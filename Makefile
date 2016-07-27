CC=mpif90
CFLAGS=-g -fcheck=all -Wall -fbacktrace
#CFLAGS=-g -O0
milagro:minimilagro.f90
	$(CC) -o milagro $(CFLAGS) $<

mini:milagro_miniapp.f90
	$(CC) -o mini $(CFLAGS) $<

strip:milagro_strip.f90
	$(CC) -o strip $(CFLAGS) $<

test:test.f90
	$(CC) -o test $(CFLAGS) $<

test2:test2.f90
	$(CC) -o test2 $(CFLAGS) $<

trip:milagro_trip.f90
	$(CC) -o trip $(CFLAGS) $<

run:milagro
	mpirun -np 2 ./milagro

runm:mini
	mpirun -np 4 ./mini

runs:strip
	mpirun -np 4 ./strip

runt:test
	mpirun -np 4 ./test

runt2:test2
	mpirun -np 4 ./test2

runp:trip
	mpirun -np 4 ./trip

clean:
	rm -f milagro mini strip test test2 trip *.o
