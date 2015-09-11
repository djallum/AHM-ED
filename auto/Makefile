# define the compiler to use:
CC = gfortran 

# compiler flags:
#   -g  adds debugging information to the executable file
#   -O0 and -03 defines the level of optimization (0 is low 3 is high)
#CFLAGS = -g -O0 
CFLAGS = -O3

#CFLAGS += -fbackslash

#define any libraries to be linked to the executable:
LIBS = -llapack

main.e: lapack.o routines.o main.o
	$(CC) $(CFLAGS) lapack.o routines.o main.o $(LIBS) -o main.e

routines.o: routines.f90
	$(CC) $(CFLAGS) -c routines.f90

lapack.o: lapack.f90
	$(CC) $(CFLAGS) -c lapack.f90

main.o: main.f90 routines.o
	$(CC) $(CFLAGS) -c main.f90

clean: 
	rm -f *.o *.mod *~