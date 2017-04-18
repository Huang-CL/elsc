DEBUG=-g
CC=mpic++
CFLAGS+= -O3 -Wall $(DEBUG)

MonteCarlo: main.o interp.o random.o spectra.o redistribute.o sphere_vr-2.o
#MonteCarlo: main.o interp.o random.o spectra.o redistribute.o sphere_vr.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^

clean:
	rm -f *.o MonteCarlo
