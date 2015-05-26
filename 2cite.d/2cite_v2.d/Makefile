ground.e: random_generator.o hamiltonian.o ground.o
	gfortran hamiltonian.o random_generator.o ground.o -llapack -o ground.e

hamiltonian.o: hamiltonian.f90
	gfortran -c hamiltonian.f90

random_generator.o : random_generator.f90
	gfortran -c random_generator.f90

ground.o : ground.f90
	gfortran -c ground.f90