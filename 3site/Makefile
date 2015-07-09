# define the compiler to use:
CC = gfortran

# compiler flags:
#   -g  adds debugging information to the executable file
#   -O0 and -03 defines the level of optimization (0 is low 3 is high)
#CFLAGS = -g -O0
CFLAGS = -O3

#define any libraries to be linked to the executable:
LIBS = -llapack

#define the executable you would like to compile
MAIN = main.e

# define the source code files needed to compile the executable
SRCS = main.f90 routines.f90

# define the object files (change the endings from .f90 to .o)
OBJS = $(SRCS:.f90=.o)

all: $(MAIN)
	@echo compiled executable: $(MAIN)

main.e: $(OBJS)
	$(CC) -o $(MAIN) $(OBJS) $(LIBS)

routines.o: routines.f90
	$(CC) $(CFLAGS) routines.f90 -c

main.o: main.f90
	$(CC) $(CFLAGS) main.f90 -c

.PHONY: clean
clean:
	$(RM) *.o

