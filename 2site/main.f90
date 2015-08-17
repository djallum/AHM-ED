program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/05/18            created code 
! P. Daley              15/07/20            rewrote LDOS calculations
! P. Daley              15/07/21            added comments
! P. Daley              15/07/23            corrected phase_PES
! P. Daley              15/08/05            added comments
!
! This code solves the 2site problem with on-site interactions and hopping. The hamiltonian matrices are found by hand and entered into the program.
! The matrices are solved using lapack. The wieght of peaks for the LDOS is found by multiplying the ground state vector by matrices representing PES and 
! IPES from both cites and computed its inner product with each eigenstate. The location of the peaks is calculated by subtracting or adding (depending
! if PES or IPES) the lowest grandpotential with the grandpotential of that eigenstate. The DOS and GIPR are calculated using the two LDOS with the 
! formulas in J. Perera and R. Wortis's paper "Energy dependence of localization with interactons and disorder: The ...". The GIPR is also weighted with 
! the equation found in the paper.
!
! The fock state basis used in order is (+ is up, - is down) :
! 00,0+,+0,0-,-0,++,--,-+,+-,02,20,+2,2+,-2,2-,22

	use routines

	implicit none

	!-------------------Input Parameters---------------------------------
	integer, parameter :: npairs = 10000000  ! size of the ensemble
	real, parameter :: t = -1                 ! nearest neighbour hoping 
	real, parameter :: delta = 12.0           ! width of disorder for the site potentials 
  	real, parameter :: U = 4.0                ! on-site interactions
 	real :: mu = U/2               ! chemical potential (half filling)
 	integer, parameter :: nbins = 600         ! number of bins for energy bining to get smooth curves
  	real, parameter :: frequency_max = 12     ! maximum energy considered in energy bining
  	real, parameter :: frequency_min = -12    ! lowest energy considered in energy bining

  	!-------------------Dependent Variables------------------------------
 	integer :: pair, i, j                        ! counters for loops
  	integer :: error                    		 ! variable for error message
  	integer :: location(1)                       ! stores the location in the grand_potential array of the lowest energy
  	integer :: bin=0                             ! index for the bin number the peak goes in
  	real, dimension(nsites) :: E                 ! site potentials
  	real, dimension(nsites,nstates) :: PES_dn_ground, PES_up_ground   ! ground state vector after an down or up photo emmision respectively
  	real, dimension(nsites,nstates) :: IPES_dn_ground, IPES_up_ground ! ground state vector after an down or up inverse photo emmision respectively
  	real, dimension(nsites,2*nstates,2) :: LDOS     ! Local DOS (LDOS(i,j,1) is the energy of state j's contribution LDOS for site i and LDOS(i,j,2) is the weight of the contribution)
  	real :: inner_product_up, inner_product_dn   ! used in the calculation of the weight of LDOS contributions
  	real :: frequency_delta                      ! step size between different energy bins
  	real, dimension(nbins,2) :: DOS=0            ! array that stores the DOS peaks and all the frequencies of the energy bins (DOS(i,1) is frequency of bin i and DOS(i,2) is the weight)
  	real :: IPR(512)                             ! contains the Inverse Participation Ratio (IPR) contributions (non weighted)
  	real, dimension(nbins) :: GIPR_num=0, GIPR_den=0, GIPR=0  ! arrays that store the numerator and denominator and the weighted IPR
  	real :: dos_sum                              ! sums the entire DOS so that it can be normalized to 1
  	real :: half_sum                             ! sums half of the DOS to find the filling 
  	character(len=50) :: filename

  	!-----------------Preparations for the main loop----------------------

  	frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

  	call random_gen_seed()
  	call transformations()    ! calls subroutines in routines.f90. It makes the lookup tables for PES and IPES.

  	call make_filename(filename,t,U,mu,delta)   ! calls subroutine in routines.f90 that assigns a filename based on the input parameters

 	open(unit=10,file=filename, status='replace', action='write', IOSTAT = error)   ! open the file that DOS and GIPR will be printed to
  	if (error/=0) then                                                              ! check if there was error when opening the file
    	write(*,*) "Error opening output file. Error number:", error                ! if there was an error print the error message
    	stop                                                                        ! exit the program 
  	end if

  	! write informartion about the code to the top of the output file
  	write(10,*) "# 2site DOS and GIPR created by main.f90 with subroutines in routines.f90"
	write(10,*) "#"
	write(10,'(A17,F6.2,2(A5,F6.2))') "# parameters: U =",U,"t =",t,"W =",delta
	write(10,'(A9,I4,A20,I12)') "# nbins =",nbins, "ensemble size =", npairs 
	write(10,*) "#"
  	
  	!---------------------------Main loop-----------------------------------

	pairs: do pair=1,npairs

		if (npairs < 100) goto 15                     ! skips the pecentage completion loop if npairs < 100 since it would cause a segmentation fault
		if (MOD(pair,npairs/100) == 0) then
     		write(*,*) nint(real(pair)/npairs*100), "%"    ! this section will print the percentage of completion.
    	end if
    	15 continue

		! zero variables for each loop

		v_ground = 0.0
	    LDOS = 0.0
	    inner_product_up = 0.0
	    inner_product_dn = 0.0

	    call site_potentials(delta,E)         ! call subroutine to assign the site potentials
	    call hamiltonian(t,U,mu,E)            ! make and solve the hamiltonians for their eigenvectors and eigenvalues

	    !-----find the ground state grand potential------------------------

	    grand_potential_ground = minval(grand_potential)   ! find the lowest grand ensemble energy

	    !-----find the corresponding eigenvector----------------

	    location = minloc(grand_potential)          ! find the location of the lowest energy  
	    v_ground = eigenvectors(location(1),:)      ! set v ground to the eigenvector corresponding to the lowest energy

	    !-----multiply v_ground by the PES and IPES matrices-------------

	    do j=1,nsites
	       	do i=1,nstates
	          	if (PES_up(j,i)==0) then
	             	PES_up_ground(j,i) = 0.0
	          	else 
	             	PES_up_ground(j,i) = v_ground(PES_up(j,i))*phase_PES_up(j,i)
	          	end if
	          	if (PES_dn(j,i)==0) then
	             	PES_dn_ground(j,i) = 0.0
	          	else 
	             	PES_dn_ground(j,i) = v_ground(PES_dn(j,i))*phase_PES_dn(j,i)
	          	end if
	          	if (IPES_up(j,i)==0) then
	             	IPES_up_ground(j,i) = 0.0
	          	else 
	             	IPES_up_ground(j,i) = v_ground(IPES_up(j,i))*phase_IPES_up(j,i)
	          	end if
	          	if (IPES_dn(j,i)==0) then
	             	IPES_dn_ground(j,i) = 0.0
	          	else 
	             	IPES_dn_ground(j,i) = v_ground(IPES_dn(j,i))*phase_IPES_dn(j,i)
	          	end if
	       	end do
	    end do 

	    ! calculate the LDOS for each site

	    do j=1,nsites
	       	do i=1,nstates
	          	inner_product_up = (dot_product(PES_up_ground(j,:),eigenvectors(i,:)))**2   ! dot product the PES up ground vector with each eigenstate
	          	inner_product_dn =  (dot_product(PES_dn_ground(j,:),eigenvectors(i,:)))**2  ! dot product the PES down ground vector with each eigenstate
	          	LDOS(j,i,1) = grand_potential_ground - grand_potential(i)                   ! location of the contribution
	          	LDOS(j,i,2) = (inner_product_up + inner_product_dn)*0.5                     ! weight of the contribution (average up and down spin components)
	       	end do
	    end do

	    do j=1,nsites
	       	do i=1,nstates
	          	inner_product_up = (dot_product(IPES_up_ground(j,:),eigenvectors(i,:)))**2
	          	inner_product_dn =  (dot_product(IPES_dn_ground(j,:),eigenvectors(i,:)))**2
	          	LDOS(j,i+nstates,1) = grand_potential(i) - grand_potential_ground           ! location of the contribution
	          	LDOS(j,i+nstates,2) = (inner_product_up + inner_product_dn)*0.5           ! weight of the contribution (average up and down spin components)
	       	end do
	    end do

	    do i=1,2*nstates
	       	bin = floor(LDOS(1,i,1)/frequency_delta) + nbins/2  + 1              ! find the bin number for energy bining
	       	DOS(bin,2) = DOS(bin,2) + (SUM(LDOS(:,i,2)))/real(nsites)            ! make the contribution to that bin
	       	if (SUM(LDOS(:,i,2)) /= 0) then
	          	IPR(i) = SUM(LDOS(:,i,2)**2)/(SUM(LDOS(:,i,2))**2)
	          	GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(SUM(LDOS(:,i,2)))/real(nsites)  ! numerator of the weighted GIPR
	          	GIPR_den(bin) = GIPR_den(bin) + (SUM(LDOS(:,i,2)))/real(nsites)         ! denominator of the weighted GIPR
	       	end if
	    end do

	end do pairs

	!-----------------Normalize DOS and print to output file--------------------------

  	dos_sum = DOS(1,2)                  ! set the sum equal to first term since the loop that sums the DOS starts at index 2.
  	DOS(1,1) = frequency_min            ! set the first DOS energy to the mininum frequency

  	do i=2,nbins                                    ! calculate sum to normalize the area under DOS to 1
     	DOS(i,1) = DOS(i-1,1) + frequency_delta     ! step the energy by frequency_delta
     	dos_sum = dos_sum + DOS(i,2)
  	end do

  	half_sum = DOS(1,2)                  ! set the sum equal to first term since the loop that sums the DOS starts at index 2.
  	do i=2,nbins/2                       ! calculate half sum
  	   half_sum = half_sum + DOS(i,2)
  	end do

  	! print additional information to the file
  	write(10,*) "# Filling:", half_sum/dos_sum
  	write(10,*) "# Filling Error:", DOS(nbins/2+1,2)/dos_sum
  	write(10,*) "#"
  	write(10,*) "#  Frequency           DOS            GIPR" 

  	do i=1,nbins
    	GIPR(i) = GIPR_num(i)/GIPR_den(i)
    	write(10,*), DOS(i,1), DOS(i,2)/dos_sum/frequency_delta, GIPR(i)
  	end do

  	close(10)

end program main
