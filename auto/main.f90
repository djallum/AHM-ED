program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/07/18            created code 
! P. Daley              15/08/17            changed filenaming
! P. Daley              15/08/18            added comments
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
	integer, parameter :: npairs=10    ! size of the ensemble
	real, parameter :: t = -1.0            ! nearest neighbour hoping 
	real, parameter :: U = 4             ! the on site interactions
	real, parameter :: delta=12.0          ! the width of the disorder
  	real, parameter :: mu = U/2            ! the chemical potential (U/2 is half filling)
  	integer, parameter :: nbins = 300             ! number of bins for energy bining to get smooth curves
	real, parameter :: frequency_max = 12         ! maximum energy considered in energy bining
	real, parameter :: frequency_min = -12        ! lowest energy considered in energy bining

  	!-------------------Dependent Variables------------------------------
  	real :: E(nsites)                   ! array to hold the site potentials
  	real, dimension(nsites,total_states) :: PESdn_MBG=0.0, PESup_MBG=0.0, IPESdn_MBG=0.0, IPESup_MBG=0.0
  	real, dimension(total_states) :: MBGvec
	real, dimension(nsites,2*total_states,2) :: LDOS=0.0
	real :: inner_prod_up=0.0, inner_prod_dn=0.0
	real :: IPR(2*total_states)=0.0
	real :: frequency_delta=0.0                        ! step size between different energy bins
	integer :: bin=0                                   ! index for the bin number the peak goes in
	real, dimension(nbins,2) :: DOS=0.0                                ! array that stores the DOS peaks and all the frequencies of the energy bins
	real, dimension(nbins) :: GIPR_num=0.0, GIPR_den=0.0, GIPR=0.0     ! arrays that store the numerator and denominator and the GIPR
	real :: dos_sum=0.0, half_sum=0.0
	integer :: i,j, pair
	integer :: g_up,g_dn,low,high,n_up,n_dn, min_up,min_dn, max_up, max_dn
	integer :: error=0                     ! variable for error message when opening file
 	integer :: location(1)=0, groundloc(2)               ! stores the location in the grand_potential array of the lowest energy 
 	character(len=50) :: filename 
 	character(len=50) :: str_npairs, str_nbins, str_ver   
 	integer :: version
 	integer :: file_count          

	frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

	call make_filename(filename,t,U,mu,delta)

	file_count = 0
	15 continue
	file_count = file_count + 1
	open(unit=10,file=filename, status='new', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be printed to
  	if (error/=0) then
    	if (filename(LEN_TRIM(filename) - 5:LEN_TRIM(filename) - 5) == '.') then
    		write(filename,'(A)') trim(adjustl(filename(1:LEN_TRIM(filename) - 4))) // "_1.dat"
    	else
    		read(filename(LEN_TRIM(filename) - 4:LEN_TRIM(filename) - 4),'(I1)') version
    		version = version + 1
    		write(str_ver,'(I1)') version
    		write(filename,'(A)') trim(adjustl(filename(1:LEN_TRIM(filename) - 5))) // trim(adjustl(str_ver)) // ".dat"
    	end if
    	if (file_count > 10) then
    		write(*,*) 'error opening output file. Error number:', error
    		stop
		end if
    	go to 15
  	end if

  	write(str_npairs,'(I15)') npairs
  	write(str_nbins, '(I10)') nbins
  	! write informartion about the code to the top of the output file
  	write(10,*) "# created by main.f90 with subroutines in routines.f90"
  	write(10,*) "#"
	write(10,'(A17,F6.2,2(A5,F6.2))') "# parameters: U=",U,"t=",t,"W=",delta
  	write(10,'(A)',advance='no') " # nbins= " // trim(adjustl(str_nbins)) 
  	write(10,*) " ensemble size=" // adjustl(trim(str_npairs)) 
  	write(10,*) "#" 

	call num_sites()
	call random_gen_seed()
	call transformations()
	call make_neighbours()
	call make_hamiltonian2(t)

	pairs: do pair=1,npairs

		if (npairs < 100) then
			write(*,*) pair
     	else 
     		if (MOD(pair,npairs/100) == 0) then
     			write(*,'(I3,A1)') nint(real(pair)/npairs*100), "%"
     		end if
    	end if

		MBGvec=0.0
    	grand_potential_ground = 0.0
    	LDOS = 0.0
    	PESdn_MBG=0.0; PESup_MBG=0.0; IPESdn_MBG=0.0; IPESup_MBG=0.0
    	grand_potential = 0
    	eigenvectors = 0
    	e_ground = 0

		call site_potentials(delta,E)
		call solve_hamiltonian1(E,U,mu)

		!-----find ground state energy------------------------

	    grand_potential_ground = minval(e_ground)   ! find the lowest grand ensemble energy

	    !-----find the corresponding eigenvector---------------

	    groundloc = minloc(e_ground)
	    g_up = groundloc(1) - 1 
	    g_dn = groundloc(2) - 1

	    location = mblock(g_up,g_dn)          ! find the location of the lowest energy  

	   	min_up = MAX(0,g_up-1)
	    max_up = MIN(nsites,g_up+1)
	    min_dn = MAX(0,g_dn-1)
	    max_dn = MIN(nsites,g_dn+1)

	    !write(*,*) g_up, g_dn

	    call solve_hamiltonian2(E,U,mu,g_up,g_dn)
	    if (g_up /= 0) call solve_hamiltonian2(E,U,mu,g_up-1,g_dn)
	    if (g_dn /= 0) call solve_hamiltonian2(E,U,mu,g_up,g_dn-1)
	    if (g_up /= nsites) call solve_hamiltonian2(E,U,mu,g_up+1,g_dn)
	    if (g_dn /= nsites) call solve_hamiltonian2(E,U,mu,g_up,g_dn+1)

	   	high = mblock(g_up,g_dn) + msize(g_up,g_dn) - 1
	   	MBGvec(mblock(g_up,g_dn):high) = eigenvectors(location(1),1:msize(g_up,g_dn))      ! set v ground to the eigenvector corresponding to the lowest energy

	    !multiply ground state vector by the matrices

	    do j=1,nsites
	       do i=1,total_states
	          if (PES_up(j,i)==0) then
	             PESup_MBG(j,i) = 0.0
	          else 
	             PESup_MBG(j,i) = MBGvec(PES_up(j,i))*phase_PES_up(j,i)
	          end if
	          if (PES_down(j,i)==0) then
	             PESdn_MBG(j,i) = 0.0
	          else 
	             PESdn_MBG(j,i) = MBGvec(PES_down(j,i))*phase_PES_down(j,i)
	          end if
	          if (IPES_up(j,i)==0) then
	             IPESup_MBG(j,i) = 0.0
	          else 
	             IPESup_MBG(j,i) = MBGvec(IPES_up(j,i))*phase_IPES_up(j,i)
	          end if
	          if (IPES_down(j,i)==0) then
	             IPESdn_MBG(j,i) = 0.0
	          else 
	             IPESdn_MBG(j,i) = MBGvec(IPES_down(j,i))*phase_IPES_down(j,i)
	          end if
	       end do
	    end do 

	    ! calculate the LDOS for all the cites
	    do n_up=min_up,max_up
			do n_dn=min_dn,max_dn
				if (n_up == min_up .and. n_dn == min_dn .and. g_up /= 0 .and. g_dn /= 0) CYCLE
				if (n_up == max_up .and. n_dn == max_dn .and. g_up /= nsites .and. g_dn /= nsites) CYCLE
				if (n_up == max_up .and. n_dn == min_dn .and. g_up /= nsites .and. g_dn /= 0) CYCLE
				if (n_up == min_up .and. n_dn == max_dn .and. g_up /= 0 .and. g_dn /= nsites) CYCLE
				if (n_up == g_up .and. n_dn == g_dn) CYCLE
				low = mblock(n_up,n_dn)
			    high = mblock(n_up,n_dn) + msize(n_up,n_dn) - 1
			    do j=1,nsites
			       	do i=low,high
			       		inner_prod_up = 0
			       		inner_prod_dn = 0
			       		if (n_up == min_up) then
			          		inner_prod_up = (dot_product(PESup_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			          	if (n_dn == min_dn) then
			          		inner_prod_dn =  (dot_product(PESdn_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			          	LDOS(j,i,1) = grand_potential_ground - grand_potential(i)              ! location of the peak
			          	LDOS(j,i,2) = (inner_prod_up + inner_prod_dn)*0.5           ! weight of the peak (average up and down spin components)
			          	inner_prod_up = 0
			       		inner_prod_dn = 0
			          	if (n_up == max_up) then
			          		inner_prod_up = (dot_product(IPESup_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			          	if (n_dn == max_dn) then
			          		inner_prod_dn =  (dot_product(IPESdn_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			         	LDOS(j,i+total_states,1) = grand_potential(i) - grand_potential_ground           ! location of the peak
			         	LDOS(j,i+total_states,2) = (inner_prod_up + inner_prod_dn)*0.5        ! weight of the peak
			       	end do
			    end do
		    end do
	    end do

	    do i=1,2*total_states
	       bin = floor(LDOS(2,i,1)/frequency_delta) + nbins/2  +1              !find the bin number for energy bining
	       DOS(bin,2) = DOS(bin,2) + (SUM(LDOS(:,i,2)))/real(nsites)
	        if (SUM(LDOS(:,i,2)) /= 0) then
	          IPR(i) = SUM(LDOS(:,i,2)**2)/(SUM(LDOS(:,i,2))**2)
	          GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(SUM(LDOS(:,i,2)))/real(nsites)  ! numerator of the weighted GIPR
	          GIPR_den(bin) = GIPR_den(bin) + (SUM(LDOS(:,i,2)))/real(nsites)         ! denominator of the weighted GIPR
	       end if
	    end do

	end do pairs

	dos_sum = DOS(1,2)
  	DOS(1,1) = frequency_min

 	do i=2,nbins                                    ! calculate sum to normalize the area under DOS to 1
    	DOS(i,1) = DOS(i-1,1) + frequency_delta
    	dos_sum = dos_sum + DOS(i,2)
  	end do

  	 half_sum = DOS(1,2)
  	do i=2,nbins/2                                    ! calculate half sum
  	   half_sum = half_sum + DOS(i,2)
  	end do

  	write(10,*) "# Filling:", half_sum/dos_sum
  	write(10,*) "# Filling Error:", DOS(nbins/2+1,2)/dos_sum
  	write(10,*) "#"
  	write(10,*) "#  Frequency           DOS             GIPR"

  	do i=1,nbins
    	GIPR(i) = GIPR_num(i)/GIPR_den(i)
    	if(DOS(i,2)/dos_sum/frequency_delta < 0.00001) then
      		GIPR(i) = 1/real(nsites)
    	end if
    	write(10,*), DOS(i,1), DOS(i,2)/dos_sum/frequency_delta, GIPR(i)
  	end do

  close(10)

end program
