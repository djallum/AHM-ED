ground.e : random_generator.o parameters.o hamiltonian.o math_setup.o groundstate.o
	gfortran hamiltonian.o random_generator.o parameters.o math_setup.o groundstate.o -llapack -g -o ground.e

random_generator.o : random_generator.f90
	gfortran -g -c random_generator.f90

parameters.o: parameters.f90
	gfortran -g -c parameters.f90

math_setup.o: math_setup.f90
	gfortran -g -c math_setup.f90

hamiltonian.o : hamiltonian.f90 math_setup.f90
	gfortran math_setup.90 -c hamiltonian.f90

groundstate.o : groundstate.f90
	gfortran -g -c groundstate.f90

