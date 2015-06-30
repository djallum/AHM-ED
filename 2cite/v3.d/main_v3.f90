program main

! for one iteration of the two site problem. All the code is contained without modules for the test. (150519). 
! A few routines are being modulated. (150521)
! this version calculates the DOS and both LDOS. No GIPR and only one iteration. Issues with the IPES matrices indices were fixed
! fock state basis in order is (+ is up, - is down) :
! 00,0+,+0,0-,-0,++,--,-+,+-,02,20,+2,2+,-2,2-,22

use random_generator
use routines

implicit none

integer, parameter :: nsites = 2
real, parameter :: t=-1 ! hopping
real, parameter :: U=10 ! on-site interaction
real, parameter :: mu=1 ! chemical potential (half-filling)
real, parameter :: delta =2.5
real :: E(nsites)
real :: H0=0,H1(2,2)=0,H2(2,2)=0,H3=0,H4=0,H5(4,4)=0,H6(2,2)=0,H7(2,2)=0,H8=0  !hamiltonian sub matrices
real :: W0=0,W1(2)=0.,W2(2)=0,W3=0.,W4=0,W5(4)=0,W6(2)=0,W7(2)=0,W8=0  !matrices for eigenvectors
real :: omega(16) ! grand esemble energies (eigen_energies - mu*number_electrons)
real :: omega_ground=0 ! the minimum grand ensemble energy
real :: eigenvector(16,16)=0 ! the 16 eigenvectors
real :: v_ground(16)=0  !ground state eigenvector
integer :: location(1)  !will store the location in the omega array of the lowest energy
integer :: PES_up(16)=0, PES_down(16)=0 ! matrices that operate on eigenstate for photo emmisions (up2 mean remove up spin from cite 2 etc.)
integer :: IPES_up(16)=0,IPES_down(16)=0! matrices for inverse photo emmisions (up mean add spin up electron)
integer :: phase_up(16), phase_down(16)! matrices to hold -1 or 1 for anticomutation
real :: LDOS1(32,2), LDOS2(32,2), DOS(32,2)
real :: PES_up_ground(16), PES_down_ground(16), IPES_up_ground(16), IPES_down_ground(16)
real :: inner_product_up(16), inner_product_down(16)
integer :: i !counter
!-----(for lapack)----------------
integer :: INFO = 0
integer :: LWORK = 16
integer :: WORK(16)

call random_gen_seed() ! seed the random number generator
call site_potentials(delta,E)

!-------(making matrices and solving eigenstates and eigenvaleus)-----------------------
H0 = 0
W0 = H0
omega(1) = W0
eigenvector(1,1) = 1

!---one electrons states--------------
H1(1,1) = E(1)
H1(2,2) = E(2)
H1(2,1) = t
H1(1,2) = t

H2 = H1

call ssyev('v','u',2,H1,2,W1,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H1 matrix. Error code:', INFO
   stop
end if

omega(2) = W1(1) - mu*1
omega(3) = W1(2) - mu*1 
eigenvector(2,2:3) = H1(1:2,1)
eigenvector(3,2:3) = H1(1:2,2)

call ssyev('v','u',2,H2,2,W2,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H2 matrix. Error code:', INFO
   stop
end if

omega(4) = W2(1) - mu*1
omega(5) = W2(2) - mu*1
eigenvector(4,4:5) = H2(1:2,1)
eigenvector(5,4:5) = H2(1:2,2)

!---two electrons states--------------

H3 = E(1) + E(2)
H4 = E(1) + E(2)
W3 = H3
W4 = H4

H5(1,1) = E(1) + E(2)
H5(2,2) = E(1) + E(2)
H5(3,3) = 2*E(1) + U
H5(4,4) = 2*E(2) + U
H5(1,3) = t
H5(1,4) = t
H5(2,3) = -t
H5(2,4) = -t
H5(3,1) = t
H5(3,2) = -t
H5(4,1) = t
H5(4,2) = -t


call ssyev('v','u',4,H5,4,W5,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H5 matrix. Error code:', INFO
   stop
end if

omega(6) = W3 - mu*2
omega(7) = W4 - mu*2
omega(8) = W5(1) - mu*2
omega(9) = W5(2) - mu*2
omega(10) = W5(3) - mu*2
omega(11) = W5(4) - mu*2
eigenvector(6,6) = 1
eigenvector(7,7) = 1
eigenvector(8,8:11) = H5(1:4,1)
eigenvector(9,8:11) = H5(1:4,2)
eigenvector(10,8:11) = H5(1:4,3)
eigenvector(11,8:11) = H5(1:4,4)

!-----three electron states---------------

H6(1,1) = 2*E(1) + E(2) + U
H6(2,2) = E(1) + 2*E(2) + U
H6(2,1) = -t
H6(1,2) = -t

H7 = H6

call ssyev('v','u',2,H6,2,W6,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H6 matrix. Error code:', INFO
   stop
end if

omega(12) = W6(1) - mu*3 
omega(13) = W6(2) - mu*3 
eigenvector(12,12:13) = H6(1:2,1)
eigenvector(13,12:13) = H6(1:2,2)

call ssyev('v','u',2,H7,2,W7,WORK,LWORK,INFO)
if (INFO /= 0) then
   write(*,*) 'Problem with Lapack for H7 matrix. Error code:', INFO
   stop
end if

omega(14) = W7(1) - mu*3 
omega(15) = W7(2) - mu*3
eigenvector(14,14:15) = H7(1:2,1)
eigenvector(15,14:15) = H7(1:2,2)

!------four electron states-------------------------

H8 = 2*E(1) + 2*E(2) + 2*U
W8 = H8

omega(16) = W8 - mu*4
eigenvector(16,16) = 1

!-----find ground state energy------------------------

omega_ground = minval(omega)   ! find the lowest grand ensemble energy

!-----find the corresponding eigenvector----------------

location = minloc(omega)  !find the location of the lowest energy  
v_ground = eigenvector(location(1),:) !set v ground to the eigenvector corresponding to the lowest energy


!----find the LDOS for cite 1 (average the magnitude of spin up and spin down)-------
!----make the PES matrices for cite 1
PES_up(1) = 2; phase_up(1) = 1
PES_up(3) = 6; phase_up(3) = -1 
PES_up(4) = 10; phase_up(4) = -1
PES_up(5) = 8;  phase_up(5) = -1 
PES_up(7) = 14;  phase_up(7) = 1
PES_up(9) = 12; phase_up(9) = 1
PES_up(11) = 13; phase_up(11) = 1
PES_up(15) = 16; phase_up(15) = -1

PES_down(1) = 4; phase_down(1) = 1
PES_down(2) = 10; phase_down(2) = 1
PES_down(3) = 9; phase_down(3) = -1 
PES_down(5) = 7; phase_down(4) = -1
PES_down(6) = 12;  phase_down(5) = -1
PES_down(8) = 14; phase_down(8) = -1
PES_down(11) = 15; phase_down(11) = 1
PES_down(13) = 16; phase_down(13) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (PES_up(i)==0) then
      PES_up_ground(i) = 0
   else 
      PES_up_ground(i) = v_ground(PES_up(i))*phase_up(i)
   end if
end do

do i =1,16
   if (PES_down(i)==0) then
      PES_down_ground(i) = 0
   else 
      PES_down_ground(i) = v_ground(PES_down(i))*phase_down(i)
   end if
end do

!do the inner product of the ground state thats been acted on with PES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(PES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(PES_down_ground,eigenvector(i,:)))**2
   LDOS1(i,1) = omega_ground - omega(i)
   LDOS1(i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
end do

!----make the IPES matrices for cite 1
IPES_up(2) = 1; phase_up(2) = 1
IPES_up(6) = 3; phase_up(6) = -1 
IPES_up(8) = 5; phase_up(8) = -1
IPES_up(10) = 4;  phase_up(10) = -1 
IPES_up(12) = 9;  phase_up(12) = 1
IPES_up(13) = 11; phase_up(13) = 1
IPES_up(14) = 7; phase_up(14) = 1
IPES_up(16) = 15; phase_up(16) = -1           !found by hand

IPES_down(4) = 1; phase_down(4) = 1
IPES_down(7) = 5; phase_down(7) = -1
IPES_down(9) = 3; phase_down(9) = -1 
IPES_down(10) = 2; phase_down(10) = 1
IPES_down(12) = 6;  phase_down(12) = -1
IPES_down(14) = 8; phase_down(14) = -1
IPES_down(15) = 11; phase_down(15) = 1
IPES_down(16) = 13; phase_down(16) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (IPES_up(i)==0) then
      IPES_up_ground(i) = 0
   else 
      IPES_up_ground(i) = v_ground(IPES_up(i))*phase_up(i)
   end if
end do

do i =1,16
   if (IPES_down(i)==0) then
      IPES_down_ground(i) = 0
   else 
      IPES_down_ground(i) = v_ground(IPES_down(i))*phase_down(i)
   end if
end do

!do the inner product of the ground state thats been acted on with IPES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(IPES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(IPES_down_ground,eigenvector(i,:)))**2
   LDOS1(i+16,1) = omega(i) - omega_ground
   LDOS1(i+16,2) = (inner_product_up(i) + inner_product_down(i))*0.5   !average the up spin and down spin magnitudes
end do

!------find the LDOS for cite 2-------------------------------------

!----make the PES matrices for cite 2
PES_up(1) = 3; phase_up(1) = 1
PES_up(2) = 6; phase_up(2) = 1 
PES_up(4) = 9; phase_up(4) = 1
PES_up(5) = 11;  phase_up(5) = -1 
PES_up(7) = 15;  phase_up(7) = -1
PES_up(8) = 13; phase_up(8) = -1
PES_up(10) = 12; phase_up(10) = 1
PES_up(14) = 16; phase_up(14) = -1

PES_down(1) = 5; phase_down(1) = 1
PES_down(2) = 8; phase_down(2) = 1
PES_down(3) = 11; phase_down(3) = 1 
PES_down(4) = 7; phase_down(4) = 1
PES_down(6) = 13;  phase_down(6) = 1
PES_down(9) = 15; phase_down(9) = 1
PES_down(10) = 14; phase_down(10) = 1
PES_down(12) = 16; phase_down(12) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (PES_up(i)==0) then
      PES_up_ground(i) = 0
   else 
      PES_up_ground(i) = v_ground(PES_up(i))*phase_up(i)
   end if
end do

do i =1,16
   if (PES_down(i)==0) then
      PES_down_ground(i) = 0
   else 
      PES_down_ground(i) = v_ground(PES_down(i))*phase_down(i)
   end if
end do

!do the inner product of the ground state thats been acted on with PES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(PES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(PES_down_ground,eigenvector(i,:)))**2
   LDOS2(i,1) = omega_ground - omega(i)
   LDOS2(i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
end do

!----make the IPES matrices for cite 2
IPES_up(3) = 1; phase_up(3) = 1
IPES_up(6) = 2; phase_up(6) = 1 
IPES_up(9) = 4; phase_up(9) = 1
IPES_up(11) = 5;  phase_up(11) = -1 
IPES_up(12) = 10;  phase_up(12) = 1
IPES_up(13) = 8; phase_up(13) = -1
IPES_up(15) = 7; phase_up(15) = -1
IPES_up(16) = 14; phase_up(16) = -1           !found by hand

IPES_down(5) = 1; phase_down(5) = 1
IPES_down(7) = 4; phase_down(7) = 1
IPES_down(8) = 2; phase_down(8) = 1 
IPES_down(11) = 3; phase_down(11) = 1
IPES_down(13) = 6;  phase_down(13) = 1
IPES_down(14) = 10; phase_down(14) = 1
IPES_down(15) = 9; phase_down(15) = 1
IPES_down(16) = 12; phase_down(16) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (IPES_up(i)==0) then
      IPES_up_ground(i) = 0
   else 
      IPES_up_ground(i) = v_ground(IPES_up(i))*phase_up(i)
   end if
end do

do i =1,16
   if (IPES_down(i)==0) then
      IPES_down_ground(i) = 0
   else 
      IPES_down_ground(i) = v_ground(IPES_down(i))*phase_down(i)
   end if
end do

!do the inner product of the ground state thats been acted on with IPES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(IPES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(IPES_down_ground,eigenvector(i,:)))**2
   LDOS2(i+16,1) = omega(i) - omega_ground                             !calculate the frequency of the peak
   LDOS2(i+16,2) = (inner_product_up(i) + inner_product_down(i))*0.5   !average the up spin and down spin magnitudes
end do


do i=1,32
   DOS(i,1) = LDOS1(i,1)
   DOS(i,2) = (LDOS1(i,2) + LDOS2(i,2))*0.5
   write(*,*), DOS(i,:)
end do


end program main
