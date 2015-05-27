program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/05/26            created code 
!
!
use routines

implicit none

integer, parameter :: npairs = 10000
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
real :: IPR(128)
integer, parameter :: nbins = 200                 ! number of bins for energy bining to get smooth curves
real, parameter :: frequency_max = 15             ! maximum energy considered in energy bining
real, parameter :: frequency_min = -15            ! lowest energy considered in energy bining
real :: frequency_delta                           ! step size between different energy bins
integer :: bin                                    ! index for the bin number the peak goes in
real, dimension(nbins,2) :: DOS                        ! array that stores the DOS peaks and all the energy bin frequencies 
real, dimension(nbins) :: GIPR_num, GIPR_den, GIPR     ! arrays that store the numerator and denominator and the GIPR
real :: sum

frequency_delta = (frequency_max - frequency_min)/nbins

call random_gen_seed()
call transformations()

GIPR = 0
GIPR_num = 0
GIPR_den = 0
DOS = 0

open(unit=10,file='3citedata1.dat', status='replace', action='write',IOSTAT = error) ! open the file that output will be printed to
if (error/=0) then
   write(*,*) 'error opening output file. Error number:', error
end if

pairs: do pair=1,npairs

eigenvectors = 0
E = 0
omega_ground = 0
omega = 0
LDOS = 0
IPR = 0

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
      inner_product_up(i) = (dot_product(PES_up_ground(j,:),eigenvectors(i,:)))**2
      inner_product_down(i) =  (dot_product(PES_down_ground(j,:),eigenvectors(i,:)))**2
      LDOS(j,i,1) = omega_ground - omega(i)    ! location of the peak
      LDOS(j,i,2) = (inner_product_up(i) + inner_product_down(i))*0.5 ! wieght of the peak
   end do
end do

do j=1,3
   do i=1,64
      inner_product_up(i) = (dot_product(IPES_up_ground(j,:),eigenvectors(i,:)))**2
      inner_product_down(i) =  (dot_product(IPES_down_ground(j,:),eigenvectors(i,:)))**2
      LDOS(j,i+64,1) = omega(i) - omega_ground     !location of the peak
      LDOS(j,i+64,2) = (inner_product_up(i) + inner_product_down(i))*0.5 !weight of the peak
   end do
end do

do i=1,128
   bin = floor(LDOS(1,i,1)/frequency_delta) + nbins/2                  !find the bin number for energy bining
   DOS(bin,2) = DOS(bin,2) + (LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))/3
   if ((LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2)) /= 0) then
      IPR(i) = (LDOS(1,i,2)**2 + LDOS(2,i,2)**2 + LDOS(3,i,2)**2)/((LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))**2) !calculated for each transition but
   end if !                                                                                                       not weighted with DOS
   GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))/3 !numerator of the wieghted GIPR
   GIPR_den(bin) = GIPR_den(bin) + (LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))/3  !denominator of the wieghted GIPR
end do

end do pairs

sum = DOS(1,2)
DOS(1,1) = frequency_min
do i=2,nbins
   DOS(i,1) = DOS(i-1,1) + frequency_delta
   sum = sum + DOS(i,2)
end do

do i=1,nbins
   GIPR(i) = GIPR_num(i)/GIPR_den(i)
   write(10,*), DOS(i,1), DOS(i,2)/sum/frequency_delta
end do

close(10)

end program main
