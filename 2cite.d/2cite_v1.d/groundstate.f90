program main

use random_gen
use parameters
use hamiltonian

implicit none


integer :: nsites=2, i=1
real, parameter :: t=-1.0 ! hopping
real, parameter :: U=100 ! on-site interaction
real, parameter :: mu=0 ! chemical potential (half-filling) 
real, parameter :: delta=1 ! width of disorder potential
! parameters for programming
integer, parameter :: npairs=5 ! number of pairs
real, parameter :: pi=3.141592653589793
real :: W1(2,2),W2(2,2),W3(6,6),W4(2,2),W5(2,2)

real, dimension(npairs,2) :: E
integer :: nstates = 4**2 


call gen_seed()
call site_potentials(npairs,nsites,E,delta)
call make_fockstates(nsites,fockstates)
call make_hamiltonian(nsites,t,U,fockstates)

write(*,*), H01

end program main


