program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/05/26            created code 
! P. Daley              15/07/13            added comments
!
! This code finds the density of states (DOS) and the general inverse participation (GIPR) for a three site system for half filling. The program makes the hamiltonian matrices for each combination of up and down electrons and diagonalizes them using LAPACK (linear algebra library). It compares the eigenenergies to find the lowest grand potential (the ground state). 
! The program using a pre made lookup table (found in routines.f90) to find the states that the ground state vector can photo emmit or inverse photo emmit too. These tables are called IPES and PES (inverse photo emmision and photo emmision spectrum). 
! It makes contributions to the DOS at the frequency (energy) equal to the difference of that state to the ground grand potential. This calculates the local density of states (LDOS). The LDOS for each of the sites is averaged and added to the total DOS. The DOS is made smooth by using energy (frequency) bins. 
! The DOS and GIPR are outputed to a data file. All the inputs must be made by changing the values of the variables as the are declared.

   use routines

   implicit none

   integer, parameter :: npairs = 10000  ! the size of the ensemble
   real, parameter :: t = 0              ! hopping term
   real, parameter :: delta = 12         ! width of disorder for the site potentials 
   real, parameter :: U = 2              ! on-site interactions
   real, parameter :: mu = U/2           ! chemical potential (half filling) 
   real :: E(3)=0                        ! site potentials
   integer :: pair=0,i=0,j=0, k=0        ! counters
   integer :: error=0                    ! variable for error message
   integer :: location(1)=0              ! stores the location in the grand_potential array of the lowest energy 
   real, dimension(3,64) :: PES_down_ground=0, PES_up_ground=0, IPES_down_ground=0, IPES_up_ground=0 ! ground state vector after being multiplied by PES or IPES matrix. PES_down_ground(x,:) is the ground state with down electron removed from site "x".
   real, dimension(3,128,2) :: LDOS=0                 ! the local density of states for each site and state for up and down electrons
   real :: inner_product_up=0, inner_product_down=0   ! containes the inner product for up and down electrons. Used to calculate the weight for DOS
   real :: IPR(128)=0                                 ! matrix that holds the temporary inverse participation ratio for each pair
   integer, parameter :: nbins = 200                  ! number of bins for energy bining to get smooth curves
   real, parameter :: frequency_max = 10              ! maximum energy considered in energy bining
   real, parameter :: frequency_min = -10             ! lowest energy considered in energy bining
   real :: frequency_delta=0                          ! step size between different energy bins
   integer :: bin=0                                   ! index for the bin number the peak goes in
   real, dimension(nbins,2) :: DOS=0                            ! array that stores the DOS peaks and all the energy bin frequencies 
   real, dimension(nbins) :: GIPR_num=0, GIPR_den=0, GIPR=0     ! arrays that store the numerator and denominator and the GIPR
   real :: sum=0                                                ! the area underneath the DOS. Used to normalize the DOS to 1

   frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

   call random_gen_seed()  ! seed the random generator (subroutines found in routines.f90)
   call transformations()  ! make the look up tables for the PES and IPES

   open(unit=10,file='3citedata1.dat', status='replace', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be written to
   if (error/=0) then
      write(*,*) 'error opening output file. Error number:', error  ! print an error message to the terminal if file doesn't open properly
   end if

   pairs: do pair=1,npairs  ! loop over each pair in the ensemble

      ! zero all the variables for each loop

      eigenvectors = 0
      E = 0
      Grand_potential_ground = 0
      Grand_potential = 0
      LDOS = 0
      IPR = 0

      call site_potentials(delta,E) ! randomly assigns the site potentials
      call hamiltonian(E,t,U,mu)    ! makes the hamiltonian matrices and solves them for the grand potentials and the different eigenvectors

      !-----find ground state energy------------------------

      Grand_potential_ground = minval(Grand_potential)   ! find the lowest grand ensemble energy

      !-----find the corresponding eigenvector----------------

      location = minloc(Grand_potential)          ! find the location of the lowest energy  
      v_ground = eigenvectors(location(1),:)      ! set v ground to the eigenvector corresponding to the lowest energy

      !multiply ground state vector by the the PES and IPES matrices
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

      ! calculate the LDOS for all the cites
      do j=1,3
         do i=1,64
            inner_product_up = (dot_product(PES_up_ground(j,:),eigenvectors(i,:)))**2
            inner_product_down =  (dot_product(PES_down_ground(j,:),eigenvectors(i,:)))**2
            LDOS(j,i,1) = Grand_potential_ground - Grand_potential(i)           ! location of the peak
            LDOS(j,i,2) = (inner_product_up + inner_product_down)*0.5           ! weight of the peak (average up and down spin components)
         end do
      end do

      do j=1,3
         do i=1,64
            inner_product_up = (dot_product(IPES_up_ground(j,:),eigenvectors(i,:)))**2   ! the inner product of the up IPES vector with a one of the eigenvectors
            inner_product_down =  (dot_product(IPES_down_ground(j,:),eigenvectors(i,:)))**2 ! the inner product of the down IPES vector with a one of the eigenvectors
            LDOS(j,i+64,1) = Grand_potential(i) - Grand_potential_ground        ! location of the peak
            LDOS(j,i+64,2) = (inner_product_up + inner_product_down)*0.5        ! weight of the peak
         end do
      end do

      do i=1,128
         bin = floor(LDOS(2,i,1)/frequency_delta) + nbins/2                !find the bin number for energy bining
         DOS(bin,2) = DOS(bin,2) + (LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))/3.0  ! make the contribution to the DOS by averaging the LDOS
         if ((LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2)) /= 0) then
            IPR(i) = (LDOS(1,i,2)**2 + LDOS(2,i,2)**2 + LDOS(3,i,2)**2)/((LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))**2) !for each transition not weighted with DOS
            GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))/3.0  ! numerator of the weighted GIPR
            GIPR_den(bin) = GIPR_den(bin) + (LDOS(1,i,2) + LDOS(2,i,2) + LDOS(3,i,2))/3.0         ! denominator of the weighted GIPR
         end if
      end do

   end do pairs

   sum = DOS(1,2)
   DOS(1,1) = frequency_min            ! set the energy of the first bin equal to the minimum energy (frequency)

   do i=2,nbins                                    ! calculate sum to normalize the area under DOS to 1
      DOS(i,1) = DOS(i-1,1) + frequency_delta
      sum = sum + DOS(i,2)
   end do

   do i=1,nbins                                    ! print the DOS and GIPR to the output file and normalize DOS at same time
      GIPR(i) = GIPR_num(i)/GIPR_den(i)
      write(10,*), DOS(i,1), DOS(i,2)/sum/frequency_delta, GIPR(i)
   end do

   close(10)

end program main
