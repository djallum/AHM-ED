program main

use routines

implicit none

integer, parameter :: npairs=10000
real, parameter :: t = 0             ! hopping term
real, parameter :: delta = 12         ! width of disorder for the site potentials 
real, parameter :: U = 0             ! on-site interactions
real, parameter :: mu = U/2          ! chemical potential (half filling) 
real :: E(4)=0                       ! site potentials
integer :: pair=0,i=0,j=0, k=0       ! counters
integer :: error=0                   ! variable for error message
integer :: location(1)=0             ! stores the location in the grand_potential array of the lowest energy 
real, dimension(4,256) :: PES_down_ground=0, PES_up_ground=0, IPES_down_ground=0, IPES_up_ground=0
real, dimension(4,512,2) :: LDOS=0
real :: inner_product_up=0, inner_product_down=0
real :: IPR(512)=0
integer, parameter :: nbins = 200                  ! number of bins for energy bining to get smooth curves
real, parameter :: frequency_max = 30               ! maximum energy considered in energy bining
real, parameter :: frequency_min = -30              ! lowest energy considered in energy bining
real :: frequency_delta=0                          ! step size between different energy bins
integer :: bin=0                                   ! index for the bin number the peak goes in
real, dimension(nbins,2) :: DOS=0                            ! array that stores the DOS peaks and all the energy bin frequencies 
real, dimension(nbins) :: GIPR_num=0, GIPR_den=0, GIPR=0     ! arrays that store the numerator and denominator and the GIPR
real :: sum=0

frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

call random_gen_seed()
call transformations()

open(unit=10,file='4citedata1.dat', status='replace', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be printed to
if (error/=0) then
   write(*,*) 'error opening output file. Error number:', error
end if

pairs: do pair=1,npairs

eigenvectors = 0
E = 0
Grand_potential_ground = 0
Grand_potential = 0
LDOS = 0
IPR = 0

call site_potentials(delta,E)
call hamiltonian(E,t,U,mu)

!-----find ground state energy------------------------

Grand_potential_ground = minval(Grand_potential)   ! find the lowest grand ensemble energy

!-----find the corresponding eigenvector----------------

location = minloc(Grand_potential)          ! find the location of the lowest energy  
v_ground = eigenvectors(location(1),:)      ! set v ground to the eigenvector corresponding to the lowest energy

!multiply ground state vector by the matrices
do j=1,4
   do i=1,256
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

! calculate the LDOS for all the cites
do j=1,4
   do i=1,256
      inner_product_up = (dot_product(PES_up_ground(j,:),eigenvectors(i,:)))**2
      inner_product_down =  (dot_product(PES_down_ground(j,:),eigenvectors(i,:)))**2
      LDOS(j,i,1) = Grand_potential_ground - Grand_potential(i)           ! location of the peak
      LDOS(j,i,2) = (inner_product_up + inner_product_down)*0.5           ! weight of the peak (average up and down spin components)
   end do
end do

do j=1,4
   do i=1,256
      inner_product_up = (dot_product(IPES_up_ground(j,:),eigenvectors(i,:)))**2
      inner_product_down =  (dot_product(IPES_down_ground(j,:),eigenvectors(i,:)))**2
      LDOS(j,i+256,1) = Grand_potential(i) - Grand_potential_ground        ! location of the peak
      LDOS(j,i+256,2) = (inner_product_up + inner_product_down)*0.5        ! weight of the peak
   end do
end do

do i=1,512
   bin = floor(LDOS(2,i,1)/frequency_delta) + nbins/2                !find the bin number for energy bining
   DOS(bin,2) = DOS(bin,2) + (LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2) + LDOS(4,i,2))/4.0
   if ((LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2) + LDOS(4,i,2)) /= 0) then
      IPR(i) = (LDOS(1,i,2)**2 + LDOS(2,i,2)**2 + LDOS(3,i,2)**2 + LDOS(4,i,2)**2) 
      IPR(i) = IPR(i)/((LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2) + LDOS(4,i,2))**2)
      GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2) + LDOS(4,i,2))/4.0  ! numerator of the weighted GIPR
      GIPR_den(bin) = GIPR_den(bin) + (LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2) + LDOS(4,i,2))/4.0         ! denominator of the weighted GIPR
   end if
end do

end do pairs

sum = DOS(1,2)
DOS(1,1) = frequency_min

do i=2,nbins                                    ! calculate sum to normalize the area under DOS to 1
   DOS(i,1) = DOS(i-1,1) + frequency_delta
   sum = sum + DOS(i,2)
end do

do i=1,nbins
   GIPR(i) = GIPR_num(i)/GIPR_den(i)
   write(10,*), DOS(i,1), DOS(i,2)/sum/frequency_delta, GIPR(i)
end do

close(10)

end program main
