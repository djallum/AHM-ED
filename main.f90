program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/05/26            created code 
!
!
use routines

implicit none

integer, parameter :: npairs = 1
real, parameter :: t = 1  !hopping term
real, parameter :: delta = 10 ! width of disorder 
real, parameter :: U = 4 ! on-site interactions
real, parameter :: mu = U/2 ! chemical potential (half filling) 
real :: E(3) ! site potentials
integer :: pair,i,j ! counter
integer :: error ! variable for error message

call random_gen_seed()

open(unit=10,file='3citedata1.dat', status='replace', action='write',IOSTAT = error) ! open the file that output will be printed to
if (error/=0) then
   write(*,*) 'error opening output file. Error number:', error
end if

pairs: do pair=1,npairs

eigenvectors = 0
E = 0
omega_ground = 0
omega = 0

call site_potentials(delta,E)
call hamiltonian(E,t,U,mu)

write(*,*), omega

end do pairs

close(10)

end program main
