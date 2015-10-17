program main

! Programmer         Date(yy/mm/dd)         Modification
! -----------        --------------         ------------
! P. Daley              15/07/18            created code 
! P. Daley              15/08/17            changed filenaming
! P. Daley              15/08/18            added comments

	use routines

	implicit none

	!-------------------Input Parameters---------------------------------
	integer, parameter :: npairs=1   	 ! size of the ensemble
	real, parameter :: t = 0.0            ! nearest neighbour hoping 
	real, parameter :: U = 0             ! the on site interactions
	real, parameter :: delta=12.0          ! the width of the disorder
  	real, parameter :: mu = U/2            ! the chemical potential (U/2 is half filling)
  	integer, parameter :: nbins = 240             ! number of bins for energy bining to get smooth curves
	real, parameter :: frequency_max = 12         ! maximum energy considered in energy bining
	real, parameter :: frequency_min = -12        ! lowest energy considered in energy bining

  	!-------------------Dependent Variables------------------------------
  	integer :: pair, i, j                        ! counters for loops
  	integer :: bin                               ! index for the bin number the peak goes in DOS
  	integer :: error                             ! variable for error message when opening file
  	integer :: n_up, n_dn                        ! counters for loops (number of up or down electrons)
  	integer :: g_up,g_dn                         ! number of electrons in the many-body ground state
  	integer :: min_up, min_dn, max_up, max_dn    ! max/min number of electrons that a many-body eigenstate that can be transitioned to can have
  	integer :: low, high                         ! range in array grand_potentials of states with a specified n_up, n_dn electrons
  	integer :: groundloc(2)                      ! stores the location in array e_ground of minimum energy (this will find g_up, g_dn)
  	integer :: location(1)                       ! stores the location in the grand_potential array of the lowest energy
  	integer :: version                           ! the version number for the data file (multiple simulations with same parameters)
 	integer :: file_count                        ! counts amount of times the program has tried to open output file (gives error if > 12)
  	real :: E(nsites)                            ! array to hold the site potentials
  	real, dimension(nstates) :: MBGvec           ! many-body ground state vector (MBG) written in fock basis (MBGvec(i) is coefficient of fock state "i")
  	real, dimension(nsites,nstates) :: PESdn_MBG, PESup_MBG   ! MBG after a down or up photo emmision (PE) respectively (PESdn_MBG(i,:) is  c_{i,dn}|Psi0> )
  	real, dimension(nsites,nstates) :: IPESdn_MBG, IPESup_MBG ! MBG after a down or up inverse photo emmision respectively (IPESdn_MBG(i,:) is  c^{dagger}_{i,dn}|Psi0> )
	real, dimension(nsites,2*nstates,2) :: LDOS      ! local density of states (LDOS(i,:) is LDOS of site i)
	real :: inner_prod_up, inner_prod_dn             ! inner products used when calculating weight of LDOS contributions (<Psi|PESdn_MBG> or <Psi|IPESdn_MBG>)
	real :: IPR(2*nstates)                               ! generalized inverse participation ratio before being ensemble averaged
	real :: frequency_delta                              ! step size between different energy bins for DOS and GIPR
	real, dimension(nbins,2) :: DOS=0                          ! array that stores the DOS (DOS(i,1) is energy of bin i and DOS(i,2) is the weight of it)
	real, dimension(nbins) :: GIPR_num=0, GIPR_den=0, GIPR=0   ! arrays that store the numerator and denominator and the ensemble averaged GIPR
	real :: dos_sum                                      ! sums the entire DOS so that area under it can be normalized to 1
	real :: half_sum                                     ! sums DOS from -W/2 to 0 to find the filling
 	character(len=50) :: filename                        ! variable that stores the output filename
 	character(len=50) :: str_npairs, str_nbins           ! strings for information at top of output file
 	character(len=50) :: str_ver                         ! string for version number of output file

  	!************************ Preparations for the main loop **********************************

	frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

	!-----------------------Assigning filename----------------------------

	call make_filename(filename,t,U,mu,delta)          ! assigns filename based on the parameters of the simulation

	file_count = 0
	15 continue
	file_count = file_count + 1                                              ! add another time it has tried to open file
	open(unit=10,file=filename, status='new', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be printed to
  	if (error/=0) then                                                       ! error means that file of that name already exists (add new version now)
    	if (filename(LEN_TRIM(filename) - 5:LEN_TRIM(filename) - 5) == '.') then                    ! check if only one version already exists
    		write(filename,'(A)') trim(adjustl(filename(1:LEN_TRIM(filename) - 4))) // "_1.dat"     ! add _1 before .dat
    	else                                                                                        ! means more then one version already exists
    		read(filename(LEN_TRIM(filename) - 4:LEN_TRIM(filename) - 4),'(I1)') version            ! find out the current version number
    		version = version + 1                                                                   ! increase the version by 1
    		write(str_ver,'(I2)') version                                                           ! convert to a string
    		write(filename,'(A)') trim(adjustl(filename(1:LEN_TRIM(filename) - 5))) // trim(adjustl(str_ver)) // ".dat" ! add new version number
    	end if
    	if (file_count > 105) then          ! if tried this many times it may be due to seperate error or already 100 versions                                  
    		write(*,*) 'error opening output file. Error number:', error  ! print error message then
    		stop
		end if
    	go to 15                                                          ! try again with updated file name
  	end if

  	!------------------Writing info to data file---------------------------

  	write(str_npairs,'(I15)') npairs           ! convert integer into a string
  	write(str_nbins, '(I10)') nbins            ! convert integer into a string

  	write(10,*) "# created by main.f90 with subroutines in routines.f90"
  	write(10,*) "#"
	write(10,'(A17,F6.2,2(A5,F6.2))') "# parameters: U=",U,"t=",t,"W=",delta
  	write(10,'(A)',advance='no') " # nbins= " // trim(adjustl(str_nbins)) 
  	write(10,*) " ensemble size=" // adjustl(trim(str_npairs)) 
  	write(10,*) "#" 

  	!------------Calling subroutines needed to start loop-------------------

	call num_sites()              ! creates the fock state basis
	call random_gen_seed()        ! seeds random number generator (used for asssigning site potentials)
	call transformations()        ! make the look up tables for c_{i,sigma} (PES) and c^{dagger}_{i,sigma} (IPES)
	call make_neighbours()        ! make lookup tables for the nearest neighbour sites of each site
	call make_hamiltonian2(t)     ! make the off diagonal terms of the hamiltonian matrices (these don't change based on site potentials)

	!****************************** Main loop *************************************

	pairs: do pair=1,npairs

		!---------Print progress of the simulation---------------
		if (npairs < 100) then                      ! if ensemble size is less then 100 it will print what pair it is on
			write(*,*) pair
     	else 
     		if (MOD(pair,npairs/100) == 0) then             
     			write(*,'(I3,A1)') nint(real(pair)/npairs*100), "%"    ! if ensemble > 100 it will print the percent of completion
     		end if
    	end if

    	!---------------Set variables to zero--------------------
		MBGvec=0.0
    	grand_potential_ground = 0.0
    	LDOS = 0.0
    	PESdn_MBG=0.0; PESup_MBG=0.0; IPESdn_MBG=0.0; IPESup_MBG=0.0
    	grand_potential = 0
    	eigenvectors = 0
    	e_ground = 0

		call site_potentials(delta,E)    ! assign the site potentials from the uniform distribution bounded between -W/2 and W/2          
		call solve_hamiltonian1(E,U,mu)  ! solve for the lowest grand potential of each hamiltonian sub-matrix

		!----------Find MBG energy------------------------

	    grand_potential_ground = minval(e_ground)   ! find the lowest grand ensemble energy

	    !-------Find n_up and n_dn of MBG-----------------

	    groundloc = minloc(e_ground)       ! find the number of up and down electrons in the lowest grand ensemble energy
	    g_up = groundloc(1) - 1            ! array indices are exactly one off from number of electrons
	    g_dn = groundloc(2) - 1

	    location = mblock(g_up,g_dn)       ! find the location of the lowest grand_potential in the array grand_potential_ground  

	    !----------Find max/min n_up,n_dn that a MBE that can be transitioned to can have (+/- one of MBG)------------------
	   	min_up = MAX(0,g_up-1)             
	    max_up = MIN(nsites,g_up+1)
	    min_dn = MAX(0,g_dn-1)
	    max_dn = MIN(nsites,g_dn+1)

	    call solve_hamiltonian2(E,U,mu,g_up,g_dn)      ! solve all eigenvectors and eigenvalues of MBG subhamiltonian

	    !---------Solve for eigenvectors and eigenvalues all sub hamiltonians with +/- 1 electron as MBG--------------------
	    if (g_up /= 0) call solve_hamiltonian2(E,U,mu,g_up-1,g_dn)
	    if (g_dn /= 0) call solve_hamiltonian2(E,U,mu,g_up,g_dn-1)
	    if (g_up /= nsites) call solve_hamiltonian2(E,U,mu,g_up+1,g_dn)
	    if (g_dn /= nsites) call solve_hamiltonian2(E,U,mu,g_up,g_dn+1)

	   	high = mblock(g_up,g_dn) + msize(g_up,g_dn) - 1                                 ! find range of indexs of fock states in the MBG's submatrix
	   	MBGvec(mblock(g_up,g_dn):high) = eigenvectors(location(1),1:msize(g_up,g_dn))   ! set MBGvec to the eigenvector corresponding to the lowest energy

	    !------------------calculate PESdn_MBG, PESup_MBG, IPESup_MBG, IPESdn_MBG------------------------------------------

	    do j=1,nsites
	       do i=1,nstates
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

	    !---------------------------calculate the LDOS for all the sites--------------------------------------------------

	    do n_up=min_up,max_up                ! loop over all possible submatrices
			do n_dn=min_dn,max_dn
				if (n_up == min_up .and. n_dn == min_dn .and. g_up /= 0 .and. g_dn /= 0) CYCLE             ! not allowed (two electrons less then MBG)
				if (n_up == max_up .and. n_dn == max_dn .and. g_up /= nsites .and. g_dn /= nsites) CYCLE   ! not allowed (two electrons more then MBG)
				if (n_up == max_up .and. n_dn == min_dn .and. g_up /= nsites .and. g_dn /= 0) CYCLE        ! not allowed (two electrons different then MBG)
				if (n_up == min_up .and. n_dn == max_dn .and. g_up /= 0 .and. g_dn /= nsites) CYCLE        ! not allowed (two electrons different then MBG)
				if (n_up == g_up .and. n_dn == g_dn) CYCLE                                                 ! not allowed (same n_up and n_dn)
				low = mblock(n_up,n_dn)                               ! lowest value of range of fock states of the submatrix
			    high = mblock(n_up,n_dn) + msize(n_up,n_dn) - 1       ! highest value of range of fock states of the submatrix
			    do j=1,nsites                  ! loop over all the sites (PESdn_MBG for c_{j,sigma} with different j's)  
			       	do i=low,high              ! loop over only the states that will be non-zero (within range of submatrix)
			       		inner_prod_up = 0
			       		inner_prod_dn = 0
			       		if (n_up == min_up) then
			          		inner_prod_up = (dot_product(PESup_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			          	if (n_dn == min_dn) then
			          		inner_prod_dn =  (dot_product(PESdn_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			          	LDOS(j,i,1) = grand_potential_ground - grand_potential(i)              ! location of the peak
			          	LDOS(j,i,2) = (inner_prod_up + inner_prod_dn)*0.5                      ! weight of the peak (average up and down spin components)
			          	inner_prod_up = 0
			       		inner_prod_dn = 0
			          	if (n_up == max_up) then
			          		inner_prod_up = (dot_product(IPESup_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			          	if (n_dn == max_dn) then
			          		inner_prod_dn =  (dot_product(IPESdn_MBG(j,low:high),eigenvectors(i,1:msize(n_up,n_dn))))**2
			          	end if
			         	LDOS(j,i+nstates,1) = grand_potential(i) - grand_potential_ground       ! location of the peak
			         	LDOS(j,i+nstates,2) = (inner_prod_up + inner_prod_dn)*0.5               ! weight of the peak (average up and down spin components)
			       	end do
			    end do
		    end do
	    end do

	    !---------------------------------Energy binning for DOS and GIPR-----------------------------------------------

	    do i=1,2*nstates
	       bin = floor(LDOS(2,i,1)/frequency_delta) + nbins/2  +1               ! find the bin number for energy bining
	       DOS(bin,2) = DOS(bin,2) + (SUM(LDOS(:,i,2)))/real(nsites)            ! add the contribution (calculated for LDOS)
	        if (SUM(LDOS(:,i,2)) /= 0) then                                     ! if LDOS is zero no GIPR contribution
	          IPR(i) = SUM(LDOS(:,i,2)**2)/(SUM(LDOS(:,i,2))**2)                ! non ensemble averaged GIPR
	          GIPR_num(bin) = GIPR_num(bin) + IPR(i)*(SUM(LDOS(:,i,2)))/real(nsites)  ! numerator of the ensemble averaged GIPR
	          GIPR_den(bin) = GIPR_den(bin) + (SUM(LDOS(:,i,2)))/real(nsites)         ! denominator of the ensemble averaged GIPR
	       end if
	    end do

	end do pairs

	!*************************** Normalize DOS and print to output file *************************************

	dos_sum = DOS(1,2)                              ! set the sum equal to first term since the loop that sums the DOS starts at index 2.
  	DOS(1,1) = frequency_min                        ! set the first DOS energy to the mininum frequency

 	do i=2,nbins                                    ! calculate sum to normalize the area under DOS to 1
    	DOS(i,1) = DOS(i-1,1) + frequency_delta     ! step the energy by frequency_delta
    	dos_sum = dos_sum + DOS(i,2)
  	end do

  	half_sum = DOS(1,2)                             ! set the sum equal to first term since the loop that sums the DOS starts at index 2.
  	do i=2,nbins/2                                    ! calculate half sum
  	   half_sum = half_sum + DOS(i,2)
  	end do

  	!------------------Print additional information to the file------------------------
  	write(10,*) "# Filling:", half_sum/dos_sum
  	write(10,*) "# Filling Error:", DOS(nbins/2+1,2)/dos_sum
  	write(10,*) "#"
  	write(10,*) "#  Frequency           DOS             GIPR"

  	!------------------Print the DOS and GIPR to the output file-----------------------
  	do i=1,nbins
    	GIPR(i) = GIPR_num(i)/GIPR_den(i)
    	if(DOS(i,2)/dos_sum/frequency_delta < 0.00001) then
      		GIPR(i) = 1/real(nsites)
    	end if
    	write(10,*), DOS(i,1), DOS(i,2)/dos_sum/frequency_delta, GIPR(i)
  	end do

  close(10)

end program
