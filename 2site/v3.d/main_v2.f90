program main

! for one iteration of the two site problem. All the code is contained without modules for the test. (150519). 
! A few routines are being modulated. (150521)
! this version calculates the DOS and both LDOS. No GIPR and only one iteration
! fock state basis in order is (+ is up, - is down) :
! 00,0+,+0,0-,-0,++,--,-+,+-,02,20,+2,2+,-2,2-,22

use random_generator

implicit none

integer, parameter :: nsites = 2
real, parameter :: t=-1 ! hopping
real, parameter :: U=5 ! on-site interaction
real, parameter :: mu=1 ! chemical potential (half-filling)
real, parameter :: delta =2
real :: E(nsites)
real :: H0=0,H1(2,2)=0,H2(2,2)=0,H3=0,H4=0,H5(4,4)=0,H6(2,2)=0,H7(2,2)=0,H8=0  !hamiltonian sub matrices
real :: W0=0,W1(2)=0.,W2(2)=0,W3=0.,W4=0,W5(4)=0,W6(2)=0,W7(2)=0,W8=0  !matrices for eigenvectors
real :: omega(16) ! grand esemble energies (eigen_energies - mu*number_electrons)
real :: omega_ground=0 ! the minimum grand ensemble energy
real :: eigenvector(16,16)=0 ! the 16 eigenvectors
real :: v_ground(16)=0  !ground state eigenvector
integer :: location(1)  !will store the location in the omega array of the lowest energy
integer :: PES_up1(16)=0,PES_up2(16)=0, PES_down1(16)=0, PES_down2(16)=0 ! matrices that operate on eigenstate for photo emmisions (up2 mean remove up spin from cite 2 etc.)
integer :: IPES_up1(16)=0,IPES_up2(16)=0,IPES_down1(16)=0, IPES_down2(16)=0 ! matrices for inverse photo emmisions (up2 mean add up spin on cite 2 etc.)
integer :: phase_up1(16),phase_up2(16), phase_down1(16), phase_down2(16) ! matrices to hold -1 or 1 for anticomutation
real :: LDOS1(32,2), LDOS2(32,2), DOS(32,2)
real :: PES_up1_ground(16), PES_down1_ground(16), IPES_up1_ground(16), IPES_down1_ground(16)
real :: PES_up2_ground(16), PES_down2_ground(16), IPES_up2_ground(16), IPES_down2_ground(16)
real :: inner_product_up(16), inner_product_down(16)
integer :: i !counter
!-----(for lapack)----------------
integer :: INFO = 0
integer :: LWORK = 16
integer :: WORK(16)
!-----(for site_potentials)-----------
real :: random1, random2

call random_gen_seed() ! seed the random number generator

!-----(site_potentials)--------------------------

call random_number(random1)
call random_number(random2)
random1 = random1 - 0.5          !centering the random numbers about 0
random2 = random2 - 0.5
E(1) = delta*max(random1,random2)
E(2) = delta*min(random1,random2)

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
PES_up1(1) = 2; phase_up1(1) = 1
PES_up1(3) = 6; phase_up1(3) = -1 
PES_up1(4) = 10; phase_up1(4) = -1
PES_up1(5) = 8;  phase_up1(5) = -1 
PES_up1(7) = 14;  phase_up1(7) = 1
PES_up1(9) = 12; phase_up1(9) = 1
PES_up1(11) = 13; phase_up1(11) = 1
PES_up1(15) = 16; phase_up1(15) = -1

PES_down1(1) = 4; phase_down1(1) = 1
PES_down1(2) = 10; phase_down1(2) = 1
PES_down1(3) = 9; phase_down1(3) = -1 
PES_down1(5) = 7; phase_down1(4) = -1
PES_down1(6) = 12;  phase_down1(5) = -1
PES_down1(8) = 14; phase_down1(8) = -1
PES_down1(11) = 15; phase_down1(11) = 1
PES_down1(13) = 16; phase_down1(13) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (PES_up1(i)==0) then
      PES_up1_ground(i) = 0
   else 
      PES_up1_ground(i) = v_ground(PES_up1(i))*phase_up1(i)
   end if
end do

do i =1,16
   if (PES_down1(i)==0) then
      PES_down1_ground(i) = 0
   else 
      PES_down1_ground(i) = v_ground(PES_down1(i))*phase_down1(i)
   end if
end do

!calculate the PES half of the LDOS for cite 1 
do i=1,16
   inner_product_up(i) = (dot_product(PES_up1_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(PES_down1_ground,eigenvector(i,:)))**2
   LDOS1(i,1) = omega_ground - omega(i)
   LDOS1(i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
end do

!----make the IPES matrices for cite 1
IPES_up1(2) = 1; phase_up1(2) = 1
IPES_up1(6) = 3; phase_up1(6) = -1 
IPES_up1(8) = 5; phase_up1(8) = -1
IPES_up1(10) = 4;  phase_up1(10) = -1 
IPES_up1(12) = 9;  phase_up1(12) = 1
IPES_up1(13) = 11; phase_up1(13) = 1
IPES_up1(14) = 7; phase_up1(14) = 1
IPES_up1(16) = 15; phase_up1(16) = -1           !found by hand

IPES_down1(4) = 1; phase_down1(4) = 1
IPES_down1(7) = 5; phase_down1(7) = -1
IPES_down1(9) = 3; phase_down1(9) = -1 
IPES_down1(10) = 2; phase_down1(10) = 1
IPES_down1(12) = 6;  phase_down1(12) = -1
IPES_down1(14) = 8; phase_down1(14) = -1
IPES_down1(15) = 11; phase_down1(15) = 1
IPES_down1(16) = 13; phase_down1(16) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (IPES_up1(i)==0) then
      IPES_up1_ground(i) = 0
   else 
      IPES_up1_ground(i) = v_ground(IPES_up1(i))*phase_up1(i)
   end if
end do

do i =1,16
   if (IPES_down1(i)==0) then
      IPES_down1_ground(i) = 0
   else 
      IPES_down1_ground(i) = v_ground(IPES_down1(i))*phase_down1(i)
   end if
end do

do i=1,16
   inner_product_up(i) = (dot_product(IPES_up1_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(IPES_down1_ground,eigenvector(i,:)))**2
   LDOS1(i+16,1) = omega(i) - omega_ground
   LDOS1(i+16,2) = (inner_product_up(i) + inner_product_down(i))*0.5   !average the up spin and down spin magnitudes
end do

!------find the LDOS for cite 2-------------------------------------

!----make the PES matrices for cite 2
PES_up2(1) = 3; phase_up2(1) = 1
PES_up2(2) = 6; phase_up2(2) = 1 
PES_up2(4) = 9; phase_up2(4) = 1
PES_up2(5) = 11;  phase_up2(5) = -1 
PES_up2(7) = 15;  phase_up2(7) = -1
PES_up2(8) = 13; phase_up2(8) = -1
PES_up2(10) = 12; phase_up2(10) = 1
PES_up2(14) = 16; phase_up2(14) = -1

PES_down2(1) = 5; phase_down2(1) = 1
PES_down2(2) = 8; phase_down2(2) = 1
PES_down2(3) = 11; phase_down2(3) = 1 
PES_down2(4) = 7; phase_down2(4) = 1
PES_down2(6) = 13;  phase_down2(6) = 1
PES_down2(9) = 15; phase_down2(9) = 1
PES_down2(10) = 14; phase_down2(10) = 1
PES_down2(12) = 16; phase_down2(12) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (PES_up2(i)==0) then
      PES_up2_ground(i) = 0
   else 
      PES_up2_ground(i) = v_ground(PES_up2(i))*phase_up2(i)
   end if
end do

do i =1,16
   if (PES_down2(i)==0) then
      PES_down2_ground(i) = 0
   else 
      PES_down2_ground(i) = v_ground(PES_down2(i))*phase_down2(i)
   end if
end do

!calculate the PES half of the LDOS for cite 1 
do i=1,16
   inner_product_up(i) = (dot_product(PES_up2_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(PES_down2_ground,eigenvector(i,:)))**2
   LDOS2(i,1) = omega_ground - omega(i)
   LDOS2(i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
end do

!----make the IPES matrices for cite 2
IPES_up2(3) = 1; phase_up2(3) = 1
IPES_up2(6) = 2; phase_up2(6) = 1 
IPES_up2(9) = 4; phase_up2(9) = 1
IPES_up2(11) = 5;  phase_up2(11) = -1 
IPES_up2(12) = 10;  phase_up2(12) = 1
IPES_up2(13) = 8; phase_up2(13) = -1
IPES_up2(15) = 7; phase_up2(15) = -1
IPES_up2(16) = 14; phase_up2(16) = -1           !found by hand

IPES_down2(5) = 1; phase_down2(5) = 1
IPES_down2(7) = 4; phase_down2(7) = 1
IPES_down2(8) = 2; phase_down2(8) = 1 
IPES_down2(11) = 3; phase_down2(11) = 1
IPES_down2(13) = 6;  phase_down2(13) = 1
IPES_down2(14) = 10; phase_down2(14) = 1
IPES_down2(15) = 9; phase_down2(15) = 1
IPES_down2(16) = 12; phase_down2(16) = 1

!multiply ground state vector by the matrices
do i =1,16
   if (IPES_up2(i)==0) then
      IPES_up2_ground(i) = 0
   else 
      IPES_up2_ground(i) = v_ground(IPES_up2(i))*phase_up2(i)
   end if
end do

do i =1,16
   if (IPES_down2(i)==0) then
      IPES_down2_ground(i) = 0
   else 
      IPES_down2_ground(i) = v_ground(IPES_down2(i))*phase_down2(i)
   end if
end do

do i=1,16
   inner_product_up(i) = (dot_product(IPES_up2_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(IPES_down2_ground,eigenvector(i,:)))**2
   LDOS2(i+16,1) = omega(i) - omega_ground
   LDOS2(i+16,2) = (inner_product_up(i) + inner_product_down(i))*0.5   !average the up spin and down spin magnitudes
end do


do i=1,32
   DOS(i,1) = LDOS1(i,1)
   DOS(i,2) = (LDOS1(i,2) + LDOS2(i,2))*0.5
   write(*,*), DOS(i,:)
end do


end program main
