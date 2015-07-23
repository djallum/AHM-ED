program main

	use routines

	implicit none

	integer, parameter :: npairs=100000
	real, parameter :: t = -1.0
	real :: E(nsites)
	real, parameter :: U=4.0
  	real, parameter :: mu = U/2
  	real, parameter :: delta=12.0
  	real, dimension(nsites,total_states) :: PES_down_ground=0.0, PES_up_ground=0.0, IPES_down_ground=0.0, IPES_up_ground=0.0
  	real, dimension(total_states) :: v_ground
	real, dimension(nsites,2*total_states,2) :: LDOS=0.0
	real :: inner_product_up=0.0, inner_product_down=0.0
	real :: IPR(2*total_states)=0.0
	integer, parameter :: nbins = 300                  ! number of bins for energy bining to get smooth curves
	real, parameter :: frequency_max = 12              ! maximum energy considered in energy bining
	real, parameter :: frequency_min = -12             ! lowest energy considered in energy bining
	real :: frequency_delta=0.0                        ! step size between different energy bins
	integer :: bin=0                                   ! index for the bin number the peak goes in
	real, dimension(nbins,2) :: DOS=0.0                                      ! array that stores the DOS peaks and all the frequencies of the energy bins
	real, dimension(nbins) :: GIPR_num=0.0, GIPR_den=0.0, GIPR=0.0     ! arrays that store the numerator and denominator and the GIPR
	real :: dos_sum=0.0, half_sum=0.0
	integer :: i,j, pair
	integer :: g_up,g_dn,low,high,n_up,n_dn
	integer :: error=0                     ! variable for error message
 	integer :: location(1)=0               ! stores the location in the grand_potential array of the lowest energy 

	frequency_delta = (frequency_max - frequency_min)/nbins   ! calculating the step size between bins for the energy bining process

	open(unit=10,file='auto.dat', status='replace', action='write',IOSTAT = error) ! open the file that DOS and GIPR will be printed to
  	if (error/=0) then
    	write(*,*) 'error opening output file. Error number:', error
  	end if

  	! write informartion about the code to the top of the output file
  	write(10,*) "# created by main.f90 with subroutines in routines.f90"
  	write(10,*) "#"
	write(10,'(A17,F6.2,2(A5,F6.2))') "# parameters: U=",U,"t=",t,"W=",delta
  	write(10,'(A9,I4,A20,I12)') "# nbins=",nbins, "ensemble size=", npairs 
  	write(10,*) "#"
  	write(10,*) "#     Frequency                     DOS                     GIPR" 

	call num_sites()
	call random_gen_seed()
	call transformations()
	call make_neighbours()
	call make_hamiltonian2(t)

	pairs: do pair=1,npairs

		if (MOD(pair,npairs/100) == 0) then
     		write(*,*) nint(real(pair)/npairs*100), "%"
    	end if

		v_ground=0.0
    	grand_potential_ground = 0.0
    	LDOS = 0.0
    	inner_product_up=0.0
    	inner_product_down=0.0
    	PES_down_ground=0.0; PES_up_ground=0.0; IPES_down_ground=0.0; IPES_up_ground=0.0
    	grand_potential = 0
    	eigenvectors = 0

		call site_potentials(delta,E)
		call solve_hamiltonian2(E,U,mu)

		!-----find ground state energy------------------------

	    grand_potential_ground = minval(grand_potential)   ! find the lowest grand ensemble energy

	    !-----find the corresponding eigenvector----------------

	    location = minloc(grand_potential)          ! find the location of the lowest energy  
	    v_ground = eigenvectors(location(1),:)      ! set v ground to the eigenvector corresponding to the lowest energy

	    !multiply ground state vector by the matrices

	    do j=1,nsites
	       do i=1,total_states
	          if (PES_up(j,i)==0) then
	             PES_up_ground(j,i) = 0.0
	          else 
	             PES_up_ground(j,i) = v_ground(PES_up(j,i))*phase_PES_up(j,i)
	          end if
	          if (PES_down(j,i)==0) then
	             PES_down_ground(j,i) = 0.0
	          else 
	             PES_down_ground(j,i) = v_ground(PES_down(j,i))*phase_PES_down(j,i)
	          end if
	          if (IPES_up(j,i)==0) then
	             IPES_up_ground(j,i) = 0.0
	          else 
	             IPES_up_ground(j,i) = v_ground(IPES_up(j,i))*phase_IPES_up(j,i)
	          end if
	          if (IPES_down(j,i)==0) then
	             IPES_down_ground(j,i) = 0.0
	          else 
	             IPES_down_ground(j,i) = v_ground(IPES_down(j,i))*phase_IPES_down(j,i)
	          end if
	       end do
	    end do 

	    ! calculate the LDOS for all the cites
	    do n_up=0,nsites
		do n_dn=0,nsites
		    do j=1,nsites
		    	low = mblock(n_up,n_dn)
		       	high = mblock(n_up,n_dn) + msize(n_up,n_dn) - 1
		       	do i=low,high
		          	inner_product_up = (dot_product(PES_up_ground(j,low:high),eigenvectors(i,low:high)))**2
		          	inner_product_down =  (dot_product(PES_down_ground(j,low:high),eigenvectors(i,low:high)))**2
		          	LDOS(j,i,1) = grand_potential_ground - grand_potential(i)              ! location of the peak
		          	LDOS(j,i,2) = (inner_product_up + inner_product_down)*0.5           ! weight of the peak (average up and down spin components)
		          	inner_product_up = (dot_product(IPES_up_ground(j,low:high),eigenvectors(i,low:high)))**2
		          	inner_product_down =  (dot_product(IPES_down_ground(j,low:high),eigenvectors(i,low:high)))**2
		         	LDOS(j,i+total_states,1) = grand_potential(i) - grand_potential_ground           ! location of the peak
		         	LDOS(j,i+total_states,2) = (inner_product_up + inner_product_down)*0.5        ! weight of the peak
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

  	do i=1,nbins
    	GIPR(i) = GIPR_num(i)/GIPR_den(i)
    	if(DOS(i,2)/dos_sum/frequency_delta < 0.00001) then
      		GIPR(i) = 1/real(nsites)
    	end if
    	write(10,*), DOS(i,1), DOS(i,2)/dos_sum/frequency_delta, GIPR(i)
  	end do

  	half_sum = DOS(1,2)
  	do i=2,nbins/2                                    ! calculate half sum
  	   half_sum = half_sum + DOS(i,2)
  	end do

  	write(*,*) "Filling:", half_sum/dos_sum

  close(10)

end program
