program main

use hamiltonian
use random_gen

implicit none

integer :: i,j
integer :: nsites=2, iteration=1
real, parameter :: t=-0 ! hopping
real, parameter :: U=100 ! on-site interaction
real, parameter :: mu=0 ! chemical potential (half-filling) 
real, parameter :: delta=1 ! width of disorder potential
! parameters for programming
integer, parameter :: npairs=1 ! number of pairs
real, parameter :: pi=3.141592653589793
real :: H00, H01(2,2),H10(2,2),H11(4,4),H20,H02,H22,H21(2,2), H12(2,2)
real, dimension(npairs,2) :: E
integer :: LWORK = 20
integer :: WORK(20)
integer :: info
real:: W1(2,2),W2(2,2),W3(4,4)

call gen_seed()
call site_potentials(nsites,E,delta)
call make_hamiltonian(E,U,t,iteration,H00,H10,H01,H20,H02,H11,H21,H12,H22)

call ssyev('v','u',2,H10,2,W1,WORK,LWORK,INFO)
call ssyev('v','u',2,H21,2,W2,WORK,LWORK,INFO)
call ssyev('v','u',4,H11,4,W3,WORK,LWORK,INFO)

write(*,*), H00,W1(1,1),W1(2,2),W1(1,1),W1(2,2),W3(1,1),W3(2,2),W3(3,3),W3(4,4),W2(1,1),W2(2,2),W2(1,1),W2(2,2),H02,H20,H22
end program main
