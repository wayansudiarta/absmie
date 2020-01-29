# Compile absmie program: mie scattering in absorbing medium

CC=gcc

absmie: complex.o nrutil.o absmie.c
	$(CC) -O -lm -o absmie.exe absmie.c complex.o nrutil.o

complex.o: complex.h complex.c
	$(CC) -c complex.c
nrutil.o: nrutil.h nrutil.c
	$(CC) -c nrutil.c

clean: 
	rm *.o 