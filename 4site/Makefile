cc = gfortran
CFLAGS = -g -O0
#CFLAGS = -O3 

main.e: routines.o main.o
	$(cc) routines.o main.o -llapack -o main.e

routines.o: routines.f90
	$(cc) $(CFLAGS) routines.f90 -c

main.o: main.f90
	$(cc) $(CFLAGS) main.f90 -c

clean: 
	rm *.o