obja  = math_setup.o
obja += random_generator.o
obja += timemachine.o
obja += parameters.o
obja += Hilbert_space_routines.o
obja += lapack.o
obja += simple_lanczos_v2.o
obja += simple_recursion.o

objg  = $(obja) Zgroundstate_driver.o
objs  = $(obja) Zspectrum_driver.o

libs  = -llapack
opts  = -O2 
#opts +=	-fbounds-check # check bounds
#opts += -fno-range-check # turn off range checking so that "complement=254" doesn't cause an overflow error
#opts +=	-g  # debug
opts +=	-pg # gprofile for timing
#opts += -fopenmp
cc = gfortran

# dependence files
%.o: %.f90
	@-rm -f $<~
	$(cc) -c $< -o $@ $(opts) $(libs)

# The final target TST
ground: $(objg)
	$(cc) $(objg) -o Zgroundstate.e $(opts) $(libs)

spectrum: $(objs)
	$(cc) $(objs) -o Zspectrum.e $(opts) $(libs)

.PHONY: clean
clean:
	@-rm -f *.o *.mod *~ *.h
