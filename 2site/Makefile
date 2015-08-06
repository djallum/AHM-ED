# define the compiler to use:
CC = gfortran

# compiler flags:
#   -g  adds debugging information to the executable file
#   -O0 and -03 defines the level of optimization (0 is low 3 is high)
#CFLAGS = -g -O0
CFLAGS = -O3

#define any libraries to be linked to the executable:
LIBS = -llapack

main.e: routines.o main.o
	$(CC) $(CFLAGS) routines.o main.o $(LIBS) -o main.e

routines.o: routines.f90
	$(CC) $(CFLAGS) -c routines.f90

main.o: main.f90
	$(CC) $(CFLAGS) -c main.f90

noint2.e: noint2.f90 routines.o
	$(CC) $(CFLAGS) noint2.f90 routines.o $(LIBS) -o noint2.e

clean: 
	rm -f *.o *.mod *~