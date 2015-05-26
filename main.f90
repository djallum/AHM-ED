program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/05/26            created code 
!
!
use routines

implicit none

real, parameter :: t = 1  !hopping term
real, parameter :: delta = 10 ! width of disorder
real :: E(3) ! site potentials
integer :: i ! counter

call random_gen_seed()
call site_potentials(delta,E)
call make_hamiltonians(E,t)

do i=1,3
   write(*,*), H01(i,:)
end do

end program main
