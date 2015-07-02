main.e: routines.o main.o
	gfortran routines.o main.o -llapack -o main.e

routines.o: routines.f90
	gfortran -c routines.f90

main.o: main.f90
	gfortran -c main.f90

clean:
	rm *.o

