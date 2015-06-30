program main

! Author: Patrick Daley          
! Date last modified : May 26th 2015
! 
! This code solves the 2cite problem with on-site interactions and hopping. The hamiltonian matrices are found by hand and entered into the program.
! The matrices are solved using lapack. The wieght of peaks for the LDOS is found by multiplying the ground state vector by matrices representing PES and 
! IPES from both cites and computed its inner product with each eigenstate. The location of the peaks is calculated by subtracting or adding (depending
! if PES or IPES) the lowest grandpotential with the grandpotential of that eigenstate. The DOS and GIPR are calculated using the two LDOS with the 
! formulas in J. Perera and R. Wortis's paper "Energy dependence of localization with interactons and disorder: The ...". The GIPR is also weighted with 
! the equation found in the paper.
!
! The fock state basis used in order is (+ is up, - is down) :
! 00,0+,+0,0-,-0,++,--,-+,+-,02,20,+2,2+,-2,2-,22

use random_generator
use routines

implicit none

integer, parameter :: nsites = 2                  ! number of sites (this must be set to two)
integer, parameter :: npairs=1000000             ! number of pairs in the ensemble
real, parameter :: t=-1                           ! hopping
real, parameter :: U=0                           ! on-site interaction
real :: mu= U/2                              ! chemical potential (half-filling)
real, parameter :: delta =6                      ! width of the disorder             
integer, parameter :: nbins = 200                 ! number of bins for energy bining to get smooth curves
real, parameter :: frequency_max = 10             ! maximum energy considered in energy bining
real, parameter :: frequency_min = -10            ! lowest energy considered in energy bining
real :: frequency_delta=0                           ! step size between different energy bins
integer :: bin=0                                    ! index for the bin number the peak goes in
real, dimension(nbins,2) :: DOS=0                        ! array that stores the DOS peaks and all the energy bin frequencies 
real, dimension(nbins) :: GIPR_num=0, GIPR_den=0, GIPR=0     ! arrays that store the numerator and denominator and the GIPR
real :: sum=0, sum_half=0                                       ! used to sum the DOS and normalize it
real :: E(nsites)=0                                 ! stores the random site potentials
real :: H0=0,H1(2,2)=0,H2(2,2)=0,H3=0,H4=0,H5(4,4)=0,H6(2,2)=0,H7(2,2)=0,H8=0    ! hamiltonian sub matrices
real :: W0=0,W1(2)=0,W2(2)=0,W3=0,W4=0,W5(4)=0,W6(2)=0,W7(2)=0,W8=0              ! matrices for eigenvalues
real :: omega(16)=0                                              ! grand potentials (eigen_energies - mu*number_electrons)
real :: omega_ground=0                                           ! the lowest grand ensemble energy
real :: eigenvector(16,16)=0                                     ! the 16 eigenvectors
real :: v_ground(16)=0                                           ! the ground state eigenvector
integer :: location(1)=0                                         ! will store the location in the omega array of the lowest energy
integer :: PES_up(2,16)=0, PES_down(2,16)=0                 ! matrices that operate on eigenstate for photo emmisions (up2 mean remove up spin from cite 2 etc.)
integer :: IPES_up(2,16)=0,IPES_down(2,16)=0                ! matrices for inverse photo emmisions (up mean add spin up electron)
integer :: phase_PES_up(2,16)=0, phase_PES_down(2,16)=0     ! matrices to hold -1 or 1 for anticomutation
integer :: phase_IPES_up(2,16)=0, phase_IPES_down(2,16)=0   ! matrices to hold -1 or 1 for anticomutation
real :: LDOS1(32,2)=0, LDOS2(32,2)=0
real :: PES_up_ground(16)=0, PES_down_ground(16)=0, IPES_up_ground(16)=0, IPES_down_ground(16)=0 ! vectors afterPES and IPES matrices have operatated on v_ground
real :: IPR(32)=0
real :: inner_product_up(16)=0, inner_product_down(16)=0
integer :: error=0             ! to store error messages
integer :: i,j,iterations !counters
real :: temp=0
!-----(for lapack)----------------
integer :: INFO = 0       ! varible that error message is stored in (zero means success)
integer :: LWORK = 16
integer :: WORK(16)

open(unit=10,file='2citedata1.dat', status='replace', action='write',IOSTAT = error)
if (error/=0) then
   write(*,*) 'error opening output file. Error number:', error
end if

call random_gen_seed() ! seed the random number generator
call transformations(PES_up,PES_down,IPES_up,IPES_down,phase_PES_up,phase_PES_down,phase_IPES_up,phase_IPES_down)
frequency_delta = (frequency_max - frequency_min)/nbins

pairs: do j=1,npairs

E = 0                   !zeroing everything for starting the next loop
eigenvector = 0
LDOS1 = 0; LDOS2 = 0; IPR = 0
omega = 0
omega_ground = 0
H0=0;H1=0;H2=0;H3=0;H4=0;H5=0;H6=0;H7=0;H8=0  
W0=0;W1=0;W2=0;W3=0;W4=0;W5=0;W6=0;W7=0;W8=0

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

!multiply ground state vector by the matrices
do i =1,16
   if (PES_up(1,i)==0) then
      PES_up_ground(i) = 0
   else 
      PES_up_ground(i) = v_ground(PES_up(1,i))*phase_PES_up(1,i)
   end if
end do

do i =1,16
   if (PES_down(1,i)==0) then
      PES_down_ground(i) = 0
   else 
      PES_down_ground(i) = v_ground(PES_down(1,i))*phase_PES_down(1,i)
   end if
end do

!do the inner product of the ground state thats been acted on with PES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(PES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(PES_down_ground,eigenvector(i,:)))**2
   LDOS1(i,1) = omega_ground - omega(i)
   LDOS1(i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
end do


!multiply ground state vector by the matrices
do i =1,16
   if (IPES_up(1,i)==0) then
      IPES_up_ground(i) = 0
   else 
      IPES_up_ground(i) = v_ground(IPES_up(1,i))*phase_IPES_up(1,i)
   end if
end do

do i =1,16
   if (IPES_down(1,i)==0) then
      IPES_down_ground(i) = 0
   else 
      IPES_down_ground(i) = v_ground(IPES_down(1,i))*phase_IPES_down(1,i)
   end if
end do

!do the inner product of the ground state thats been acted on with IPES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(IPES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(IPES_down_ground,eigenvector(i,:)))**2
   LDOS1(i+16,1) = omega(i) - omega_ground
   LDOS1(i+16,2) = (inner_product_up(i) + inner_product_down(i))*0.5   !average the up spin and down spin magnitudes
end do

!------find the LDOS for cite 2------------------------------------

!multiply ground state vector by the matrices
do i =1,16
   if (PES_up(2,i)==0) then
      PES_up_ground(i) = 0
   else 
      PES_up_ground(i) = v_ground(PES_up(2,i))*phase_PES_up(2,i)
   end if
end do

do i =1,16
   if (PES_down(2,i)==0) then
      PES_down_ground(i) = 0
   else 
      PES_down_ground(i) = v_ground(PES_down(2,i))*phase_PES_down(2,i)
   end if
end do

!do the inner product of the ground state thats been acted on with PES with each eigenstate
do i=1,16
   inner_product_up(i) = (dot_product(PES_up_ground,eigenvector(i,:)))**2
   inner_product_down(i) =  (dot_product(PES_down_ground,eigenvector(i,:)))**2
   LDOS2(i,1) = omega_ground - omega(i)
   LDOS2(i,2) = (inner_product_up(i) + inner_product_down(i))*0.5
end do

!multiply ground state vector by the matrices
do i =1,16
   if (IPES_up(2,i)==0) then
      IPES_up_ground(i) = 0
   else 
      IPES_up_ground(i) = v_ground(IPES_up(2,i))*phase_IPES_up(2,i)
   end if
end do

do i=1,16
   if (IPES_down(2,i)==0) then
      IPES_down_ground(i) = 0
   else 
      IPES_down_ground(i) = v_ground(IPES_down(2,i))*phase_IPES_down(2,i)
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
   bin = floor(LDOS2(i,1)/frequency_delta) + nbins/2 !find the bin number for energy bining
   DOS(bin,2) = DOS(bin,2) + (LDOS1(i,2) + LDOS2(i,2))*0.5
   if ((LDOS1(i,2) + LDOS2(i,2))**2 /= 0) then
      IPR(i) = (LDOS1(i,2)**2 + LDOS2(i,2)**2)/((LDOS1(i,2) + LDOS2(i,2))**2) !calculated for each transition but not weighted with DOS
   end if
   GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(LDOS1(i,2) + LDOS2(i,2))*0.5  !numerator of the wieghted GIPR
   GIPR_den(bin) = GIPR_den(bin) + (LDOS1(i,2) + LDOS2(i,2))*0.5  !denominator of the wieghted GIPR
end do
end do pairs


sum = DOS(1,2)
DOS(1,1) = frequency_min

do i=2,nbins
   DOS(i,1) = DOS(i-1,1) + frequency_delta
   sum = sum + DOS(i,2)
end do

!sum_half = DOS(1,2)
!do i=2,nbins/2
!   sum_half = sum_half + DOS(i,2) 
!end do

!write(*,*)'Filling:', sum_half/sum - DOS(nbins/2,2)/sum/2
!write(*,*) 'Error range:', DOS(nbins/2,2)/sum/2 

do i=1,nbins
   GIPR(i) = GIPR_num(i)/GIPR_den(i)
   write(10,*), DOS(i,1), DOS(i,2)/sum/frequency_delta, GIPR(i)
end do

close(10)

end program main
