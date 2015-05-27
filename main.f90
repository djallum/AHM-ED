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
integer :: location(1) ! will store the location in the omega array of the lowest energy
real, dimension(3,64) :: PES_down_ground, PES_up_ground, IPES_down_ground, IPES_up_ground
real, dimension(3,128,2) :: LDOS
real :: inner_product_up(64), inner_product_down(64)

call random_gen_seed()
call transformations()

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

!-----find ground state energy------------------------

omega_ground = minval(omega)   ! find the lowest grand ensemble energy

!-----find the corresponding eigenvector----------------

location = minloc(omega)  !find the location of the lowest energy  
v_ground = eigenvectors(location(1),:) !set v ground to the eigenvector corresponding to the lowest energy

!multiply ground state vector by the matrices
do j=1,3
   do i=1,64
      if (PES_up(j,i)==0) then
         PES_up_ground(j,i) = 0
      else 
         PES_up_ground(j,i) = v_ground(PES_up(j,i))*phase_PES_up(j,i)
      end if
       if (PES_down(j,i)==0) then
         PES_down_ground(j,i) = 0
      else 
         PES_down_ground(j,i) = v_ground(PES_down(j,i))*phase_PES_down(j,i)
      end if
      if (IPES_up(j,i)==0) then
         IPES_up_ground(j,i) = 0
      else 
         IPES_up_ground(j,i) = v_ground(IPES_up(j,i))*phase_IPES_up(j,i)
      end if
      if (IPES_down(j,i)==0) then
         IPES_down_ground(j,i) = 0
      else 
         IPES_down_ground(j,i) = v_ground(IPES_down(j,i))*phase_IPES_down(j,i)
      end if
   end do
end do 

LDOS = 0
! calculate the LDOS for all the cites
do j=1,3
   do i=1,64
      inner_product_up(i) = (dot_product(PES_up_ground(j),eigenvectors(i,:)))**2
      inner_product_down(i) =  (dot_product(PES_down_ground(j),eigenvectors(i,:)))**2
      LDOS(j,i,1) = omega_ground - omega(i)
      LDOS(j,i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
   end do
end do

do j=1,3
   do i=1,64
      inner_product_up(i) = (dot_product(IPES_up_ground(j),eigenvectors(i,:)))**2
      inner_product_down(i) =  (dot_product(IPES_down_ground(j),eigenvectors(i,:)))**2
      LDOS(j,i+64,1) = omega_ground - omega(i)
      LDOS(j,i+64,2) = (inner_product_up(i) + inner_product_down(i))*0.5
   end do
end do

end do pairs

close(10)

end program main
