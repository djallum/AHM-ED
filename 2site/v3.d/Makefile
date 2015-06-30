main.e: main.o random_generator.o modules.o
	gfortran random_generator.o modules.o main.o -llapack -o main.e

random_generator.o: random_generator.f90
	gfortran -c random_generator.f90

modules.o: modules.f90
	gfortran -c modules.f90

main.o: main.f90
	gfortran -c main.f90