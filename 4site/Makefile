# define the compiler you want to use
CC = gfortran

# compiler flags:
#CFLAGS = -g -O0
CFLAGS = -O3 

main.e: routines.o main.o
	$(CC) routines.o main.o -llapack -o main.e

routines.o: routines.f90
	$(CC) $(CFLAGS) routines.f90 -c

main.o: main.f90
	$(CC) $(CFLAGS) main.f90 -c

solver: solver.o
	$(CC) solver.o -llapack -o solver

solver.o: solver.f90
	$(CC) $(CFLAGS) solver.f90 -c

compare: compare.o
	$(CC) compare.o -o compare

compare.o: compare.f90
	$(CC) $(CFLAGS) compare.f90 -c

nonint.e: routines.o nonint.o
	$(CC) routines.o nonint.o -llapack -o nonint.e

nonint.o: nonint.f90
	$(CC) $(CFLAGS) nonint.f90 -c

trans: routines.o trans.o
	$(CC) $(CFLAGS) routines.o trans.o -llapack -o trans

trans.o: trans.f90
	$(CC) $(CFLAGS) -c trans.f90

clean: 
	rm -f *.o *~ *.mod