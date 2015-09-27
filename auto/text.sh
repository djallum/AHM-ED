#!/bin/bash

nsites=8
nstates="$(echo "$((4**$nsites))")"

echo "module routines

	!    %-------------------------------------------------------------------------------------%
	!    |  Fock states are represented in binary. Each one is represented with two integers:  |
	!    |       - First integer is assosiated with up electrons                               |
	!    |       - Second integer is assosiated with down electrons                            |
	!    |  An example: The two integers are 7, 22 and it's an 8 site simulation               |
	!    |       7 in binary:  00000111                                                        |
	!    |      22 in binary:  00010110                                                        |
	!    |  This means the fock state is (+ is up - is down and 2 is both):                    |	
	!    |                     000-022+                                                        |
	!    %-------------------------------------------------------------------------------------%

	use lapack

	implicit none

	!--------------Other Parameters and Variables-------------------------------
	integer, parameter :: nsites = $nsites             ! number of states
	integer, parameter :: nstates = $nstates           ! size of the fock state (FS) basis
	real :: grand_potential_ground=0.0                 ! the lowest grand ensemble energy
	integer, dimension(2,nstates) :: fock_states       ! array that stores each FS represented in binary (see above)
	integer, dimension(0:2**nsites) :: states_order    ! orders the up part of FS and down part of FS before combining them
	integer :: ne                                      ! number of electrons (counts electrons in various routines)
	integer :: tot_states_up                         ! total number of states with no down electrons
	integer, dimension(nsites,nstates) :: PES_down, PES_up   ! lookup tables for PES of many-body groundstate (MBG) transformations (ex. c|Psi0>)
	integer, dimension(nsites,nstates) :: IPES_down, IPES_up ! lookup tables for PES of MBG transformations (ex. c^{dagger}|Psi0>)
	integer, dimension(nsites,nstates) :: phase_PES_down, phase_PES_up    ! to get anticommutation sign right
	integer, dimension(nsites,nstates) :: phase_IPES_down, phase_IPES_up  ! to get anticommutation sign right

	!------------------Hamiltonian Submatrices---------------------------------
	"
for ((n_up=0; n_up<=nsites; n_up++)); do
	echo -n "	real :: "
	for ((n_dn=0; n_dn<=nsites; n_dn++)); do
		term1="$(echo -e "$nsites \n $n_up"| ./math.e)"
		term2="$(echo -e "$nsites \n $n_dn"| ./math.e)"
		final="$(echo "$(($term1*$term2))")"
		echo -n "H$n_up$n_dn($final,$final)"
		if [[ $n_dn != $nsites ]]; then
			echo -n ","
		fi
		if [[ "$final" -gt "$max" ]]; then
			max=$final
		fi
	done
	echo ""
done

echo "
	real, dimension(nstates) :: grand_potential          ! grand potentials (eigenenergies - mu*number electrons)
	real, dimension(0:nsites,0:nsites) :: e_ground       ! array of the lowest grand potential (Gpot) of each submatrix (e_ground(i,j) is lowest Gpot of Hij) 
	real, dimension(nstates,$max) :: eigenvectors        ! the many body eigenvectors (MBE) only coefficients of basis states with same n_up,n_dn as it
	integer, dimension(0:nsites) :: block, temp_block    ! block(i) is index of the first state with i electrons (when only looking at states with n_dn=0)
	integer, dimension(0:nsites) :: nstates_up           ! nstates_up(i) is the number of states with i up electrons and n_dn = 0   
	integer, allocatable, dimension(:,:) :: neighbours   ! neighbours(i,:) is all the site that are nearest neighbours to site i
	integer, dimension(0:nsites,0:nsites) :: msize       ! msize(i,j) is size of submatrix with n_up=i and n_dn=j
	integer, dimension(0:nsites,0:nsites) :: mblock      ! mblock(i,j) is index in fock_states array of first state with n_up=i,n_dn=j

contains

!************************************************************************************

subroutine make_filename(filename,t,U,mu,delta)

	! assigns a filename for the output file based on parameters of the simulation. Naming convention found in README in the folder data

	implicit none

	character(len=50), intent(out) :: filename                ! output data file name
	character(len=10) :: mu_str, t_str, U_str, W_str, s_str   ! string to hold values of parameters
	real, intent(in) :: t,U,mu,delta                          ! values of the parameters

	!--------------Convert the parameters to strings--------------------
	write(mu_str,'(F4.1)') mu 
	write(t_str,'(I2)') nint(t)
	write(W_str,'(I2)') nint(delta) 
	write(U_str,'(I2)') nint(U)
	write(s_str,'(I1)') nsites

	!-------------Contruct the file name from strings-------------------
	write(filename,'(A)') trim(adjustl('data/a')) // trim(adjustl(s_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // '-dos+ipr_t' // trim(adjustl(t_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // 'U' // trim(adjustl(U_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // 'W' // trim(adjustl(W_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // 'mu' // trim(adjustl(mu_str)) // '.dat'

end subroutine make_filename

!************************************************************************************

subroutine random_gen_seed()

	! seeds the randome number generator (used for assigning site potentials) with system clock

	implicit none

  	integer :: seed_size         ! the size of the random generators seed
	integer :: clock             ! the system clock (used to seed random generator) 
	integer :: i                 ! counter for loop
	integer, allocatable, dimension(:) :: seed    ! will store our seed

	call random_seed(size = seed_size)       ! find the size of the random generators seed
	allocate(seed(seed_size))                ! make our seed that same size
	call system_clock(count = clock)         ! find the system time
	seed=clock + 37*(/(i-1,i=1,seed_size)/)  ! create our seed using the clock
	call random_seed(put=seed)               ! seed the random generator with our seed

	deallocate(seed)

end subroutine random_gen_seed

!************************************************************************************

subroutine site_potentials(delta,E)

  	implicit none

  	real, intent(in) :: delta                ! width of the disorder (distribution from which site potentials are chosen)
  	real, intent(out) :: E(nsites)           ! the site potentials
  	real :: random                           ! random number used in intermediate step
  	integer :: i                             ! counter for loop

  	do i=1,nsites
     	call random_number(random)           ! gives a random number between 0 and 1
     	E(i) = delta*(random - 0.5)          ! centers the random numbers about 0
  	end do

end subroutine site_potentials
	
!************************************************************************************

subroutine num_sites()

	!    %----------------------------------------------------------------------------------------------%
	!    | Creates and orders the fock state basis. As an example the order for 2 site system is:       |
	!    |                                                                                              |
	!    |                                                 i=                                           |
	!    |                   | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |   |
	!    |   ----------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|   |
	!    |   fock_states(1,i)| 0 | 0 | 0 | 0 | 1 | 2 | 1 | 2 | 1 | 2  | 1  | 2  | 3	 | 3  | 3  | 3  |   |
	!    |   ----------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|   |
	!    |   fock_states(2,i)| 0 | 1 | 2 | 3 | 0 | 0 | 1 | 1 | 2 | 2  | 3  | 3  | 0  | 1  | 2  | 3  |   |
	!    |                                                                                              |
	!    %----------------------------------------------------------------------------------------------%
	!
	!    * details of the order are not particularily important and specifics can be found by working through source code
	!
	!    It first looks at only on integer (up electrons) and orders them by ne in increasing order. Then uses this ordered set
	!    of integers (variable is called order_states) to create the fock_states. The order in which each sub block of states with 
	!	 same number of electrons shows how it uses the list.
	!    
	!    Example for nsites=4
	!    	H00,H01,H02,H03,H04,H10,H11,H12,H13,H14,H20.......
	!
	!    It would move up by block of integers with n_up electrons (assigning them to fock_states(1,i)) then and cycle through the entire set
	!    (assigning fock_states(2,i) while keeping n_up in same block. (See example in the box above)

	implicit none

	integer :: i, j, istate, isite  ! counters for loops
	integer :: max_electrons        ! maximum amount of spin-up electrons (same as nsites)

	tot_states_up = 2**nsites     ! total number of fock states with n_dn = 0
	max_electrons = nsites 

	! enter in intial variables before starting loops
	block = 0  
    nstates_up = 0
    block(0) = 1
    nstates_up(0) = 1

    ! contruct the block array
	do ne = 1, max_electrons
   		nstates_up(ne) = choose(nsites,ne)          ! number of fock states with n_up=ne (n_dn=0). Choose function defined at bottom of module
   		block(ne) = block(ne-1) + nstates_up(ne-1)  ! start of ne block is the ne-1 block plus amount of ne-1 states there are
	end do

	temp_block = block                              ! temporary array that can be modified without destroying original block array
 	
 	! %-----------------------------------------------------------------------%
 	! |	Order intergers by number of ones (electrons) they have in binary     | 
 	! |                                                                       | 
 	! | Example for nsite=4:  0, 1,2,4,8, 3,5,6,9,10,12, 7,11,13,14, 15)      | 
 	! |                      |0||   1   ||      2      |     3     | 4 |      | 
 	! |                                                                       | 
 	! %-----------------------------------------------------------------------% 

	do istate=0,tot_states_up-1
		ne = 0
		do isite=1,nsites
			ne = ne + ibits(istate,isite-1,1)     ! count the number of one of that binary number (up electrons of the state)
		end do 
		states_order(temp_block(ne)) = istate     ! place it at the lowest index of the block with the same n_up
		temp_block(ne) = temp_block(ne) + 1       ! raise the lowest index of that block by 1 (don't overwrite what you just put in)
	end do

	! Uses the ordered list (states_order) created in last loop to make the fock_states
	
	istate = 1
	do ne=0,max_electrons
		do j=1,tot_states_up
			do i=block(ne), block(ne) + nstates_up(ne) - 1
				fock_states(1,istate) = states_order(i)
				fock_states(2,istate) = states_order(j)
				istate = istate + 1
			end do
		end do
	end do

end subroutine num_sites

!************************************************************************************

subroutine transformations()

	!  %---------------------------------------------------------------------------------------------------%
	!  |  Creates the lookup tables needed to calculate c_{i,sigma}|Psi0> and cc_{i,sigma}^{dagger}|Psi0>  |
	!  |       - Psi0 is the many-body ground state                                                        |
	!  |       - c_{i,sigma} is removable of electron of spin sigma from site i                            |
	!  |       - c_{i,sigma}^{dagger} is addition of electron of spin sigma on site i                      |
	!  |  Shows which fock state each of the basis vectors will transfer to after the specified removable  |
	!  |  or addition of up/dn electron from specified site.                                               |
	!  |                                                                                                   |
	!  |  PES_up(j,i) is the state (as in its index within the basis) that after removing an up electron   |
	!  |  from site j of it, will become state i.                                                          |
	!  |                                                                                                   |
	!  |  It does this by adding a '1' (binary) to the site then checking what fock state it now became    |
	!  |  by compairing the new integers to the FS basis. To get anticommutation sign right it counts the  |
	!  |  amount of up and dn electrons starting at the site the electron was added and going to           |
	!  |  site=nsites. This gives the number of anti-commutations that the creation operator would have    |
	!  |  to do.                                                                                           | 
	!  |                                                                                                   |
	!  |  The program then uses the PES tables to calculate the IPES tables since they are opposites       |
	!  %---------------------------------------------------------------------------------------------------%

	implicit none
	
	integer :: i, j                           ! counters for loops
	integer :: new_state(2)                   ! new fock state integers after the transition 
	integer :: new_index                      ! index of the new state in the array fock_states
	integer :: position, isite                ! both counters for loops

	!------------------Zero all the variables-------------------------
	PES_up = 0; PES_down = 0
    IPES_up = 0; IPES_down = 0
    phase_PES_up = 0; phase_PES_down = 0
    phase_IPES_up = 0; phase_IPES_down = 0

	!------------------Make the PES_up tables-------------------------

	do position = 1,nsites                                         ! loop over all the sites (make PES_up for PE from each site)
		do i=1,nstates                                             ! loop over each state
			ne = 0
			if (ibits(fock_states(1,i),position-1,1) == 1) then    ! if there is an electron on that site it can't be the result of PE
				PES_up(position,i) = 0                             ! zero everything then because it's not possible
				phase_PES_up(position,i) = 0                       ! zero everything then because it's not possible
			else
				do isite=position,nsites                           ! loop over all sites greater then site of PE (count number of anti-commutations)
					ne = ne + ibits(fock_states(1,i),isite-1,1)    ! count the number of up electrons it will have to commute with to be removed
				end do
				do isite=position,nsites
					ne = ne + ibits(fock_states(2,i),isite-1,1)    ! count the number of dn electrons it will have to commute with to be removed
				end do
				if (MOD(ne,2) == 0) then 
					phase_PES_up(position,i) = 1                   ! if it had to do even number of anti-commutations it is positive
				else
					phase_PES_up(position,i) = -1                  ! if it had to do odd number of anti-commutations it is positive
				end if
				new_state(1) = ibset(fock_states(1,i),position-1)  ! add up electron to that site 
				new_state(2) = fock_states(2,i)                    ! the down portion remains the same
				do j=1,nstates
					if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then 
						new_index = j                              ! find the index of the new state by compairing it to the entire FS basis
					end if
				end do
				PES_up(position,i) = new_index                     ! record the state that will when PE will become state i
			end if
		end do
	end do

	!------------------Make the PES_dn tables-------------------------

	do position = 1,nsites
		do i=1,nstates
			ne = 0
			if (ibits(fock_states(2,i),position-1,1) == 1) then    ! if there is an electron on that site it can't be the result of PE
				PES_down(position,i) = 0                           ! zero everything then because it's not possible
				phase_PES_down(position,i) = 0                     ! zero everything then because it's not possible
			else
				do isite=position+1,nsites                         ! +1 sicne the order is dn,up so wouldn't commute with the up electron on site=position
					ne = ne + ibits(fock_states(1,i),isite-1,1)    ! count the number of up electrons it will have to commute with to be removed
				end do
				do isite=position,nsites
					ne = ne + ibits(fock_states(2,i),isite-1,1)    ! count the number of dn electrons it will have to commute with to be removed
				end do
				if (MOD(ne,2) == 0) then 
					phase_PES_down(position,i) = 1                 ! if it had to do even number of anti-commutations it is positive
				else
					phase_PES_down(position,i) = -1
				end if
				new_state(2) = ibset(fock_states(2,i),position-1)  ! add dn electron to that site 
				new_state(1) = fock_states(1,i)                    ! the up portion remains the same
				do j=1,nstates
					if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
						new_index = j                               ! find the index of the new state by compairing it to the entire FS basis
					end if
				end do
				PES_down(position,i) = new_index                    ! record the state that will when PE will become state i
			end if
		end do
	end do

	!-------Find the IPES tables----------

	sites: do j=1,nsites  
     	do i=1,nstates
        	if (PES_down(j,i) /= 0) then
           		phase_IPES_down(j,PES_down(j,i)) = phase_PES_down(j,i)
           		IPES_down(j,PES_down(j,i)) = i
        	end if
        	if (PES_up(j,i) /= 0) then
           		IPES_up(j,PES_up(j,i)) = i
           		phase_IPES_up(j,PES_up(j,i)) = phase_PES_up(j,i)
        	end if
     	end do
 	end do sites

end subroutine transformations

!************************************************************************************

subroutine make_neighbours()

! makes matrix that containes information about which sites are nearest neighbours
! neighbours(i,:) is a list of all the neighbours of site i. Each site has 4 nearest neighbours normally.

	if (nsites == 8) then                  ! betts lattice for 8 sites
		allocate(neighbours(nsites,4))
		neighbours(1,1) = 8; neighbours(1,2) = 7; neighbours(1,3) = 5; neighbours(1,4) = 3
		neighbours(2,1) = 7; neighbours(2,2) = 8; neighbours(2,3) = 3; neighbours(2,4) = 5
		neighbours(3,1) = 1; neighbours(3,2) = 2; neighbours(3,3) = 4; neighbours(3,4) = 6
		neighbours(4,1) = 5; neighbours(4,2) = 3; neighbours(4,3) = 8; neighbours(4,4) = 7
		neighbours(5,1) = 2; neighbours(5,2) = 1; neighbours(5,3) = 6; neighbours(5,4) = 4
		neighbours(6,1) = 3; neighbours(6,2) = 5; neighbours(6,3) = 7; neighbours(6,4) = 8
		neighbours(7,1) = 4; neighbours(7,2) = 6; neighbours(7,3) = 1; neighbours(7,4) = 2
		neighbours(8,1) = 6; neighbours(8,2) = 4; neighbours(8,3) = 2; neighbours(8,4) = 1
	end if

	if (nsites == 4) then                  ! linear 4 site lattice
		allocate(neighbours(nsites,2))
		neighbours(1,1) = 4; neighbours(1,2) = 2
		neighbours(2,1) = 1; neighbours(2,2) = 3
		neighbours(3,1) = 2; neighbours(3,2) = 4
		neighbours(4,1) = 3; neighbours(4,2) = 1
	end if

	if (nsites == 2) then                 ! linear 2 site lattice
		allocate(neighbours(nsites,1))
		neighbours(1,1) = 2; neighbours(2,1) = 1
	end if

end subroutine make_neighbours

!************************************************************************************

subroutine matrix_sizes()

	!  %------------------------------------------------------------------------------%
	!  |  This subroutine makes the array matrix_sizes which contains the dimensions  |
	!  |  of all the Hamiltonian submatrices                                          | 
	!  |                                                                              |
	!  |  matrix_sizes(i,j) is the size of the matrix for FS with i up electrons      |
	!  |  and j down electrons.                                                       |
	!  %------------------------------------------------------------------------------%

	implicit none

	integer :: n_up,n_dn

	do n_up=0,nsites
		do n_dn=0,nsites
			msize(n_up,n_dn) = choose(nsites,n_up)*choose(nsites,n_dn)
			if (n_dn == 0 .and. n_up == 0) then
				mblock(n_up,n_dn) = 1
			else if (n_dn == 0) then
				mblock(n_up,n_dn) = mblock(n_up-1,nsites) + msize(n_up-1,nsites)
			else 
				mblock(n_up,n_dn) = mblock(n_up,n_dn-1) + msize(n_up,n_dn-1)
			end if
		end do
	end do

end subroutine matrix_sizes

!************************************************************************************

subroutine make_hamiltonian2(t)
		
	!  %--------------------------------------------------------------------------------%	
	!  |  This routine makes the hamiltonians off diagonal terms. On diagonal terms     | 
	!  |  are added during each loop.                                                   |
	!  |                                                                                |
	!  |  Any terms in the Hamiltonians that are off the main diagonal should be added  |
	!  |  here. Terms that don't depend on site potentials can also be added here.      |
	!  |                                                                                |
	!  |  The program loops over each submatrix and every state within them. It loops   |
	!  |  over each site and checks if there is an up electron and if so it then checks |
	!  |  which of it's neighbouring sites don't have one. It then removes the '1' from |
	!  |  that site and adds it to a neighbouring one and counts the amount of electrons|
	!  |  between the two sites (to get anticommutation). It then finds the indexs of   |
	!  |  the original state and the new state and at a +/- t to that column/row of the |
	!  |  submatrix. It then repeats the process for the down electrons.                |
	!  %--------------------------------------------------------------------------------%
		
	implicit none

	real, intent(in) :: t                    ! the hopping integral
	integer :: istate, isite, i, j, y        ! counters for loops
	integer :: inbr                          ! site number of the neighbour
	integer :: new_state(2)                  ! new state after a hopping
	integer :: phase                         ! the number that keeps track of anti-commutations (+ or -)
	integer :: ne                            ! counts the number of electrons (used to calculate the phase)
	integer :: trans_site(2)                 ! holds the site number of the starting and ending site of a hopping (needed to find number of electrons)
	integer :: new_index                     ! the column index of the new FS after the hopping
	integer :: state_index                   ! the row index of the old FS before the hopping
	integer :: n_up,n_dn                     ! the number of up ad down electrons (used to loop over each submatrix)

	!------------------Zero the matrices------------------------
	"

for ((n_up=0;n_up<=nsites;n_up++)); do
	echo -n "	"
	for ((n_dn=0;n_dn<=nsites;n_dn++)); do
		echo -n "H$n_up$n_dn=0.0"
		if [[ n_dn != $nsites ]]; then
			echo -n "; "
		fi
	done
	echo ""
done

echo "
	call matrix_sizes()            ! contruct the array that tells the sizes of the Hamiltonian submatrices

	do n_up = 0,nsites             ! loop over all submatrices
	do n_dn = 0,nsites             ! loop over all submatrices
	do istate = mblock(n_up,n_dn),mblock(n_up,n_dn) + msize(n_up,n_dn)-1  ! loop over all the states in each submatrix
		do isite = 1,nsites                                               ! loop over each site of the state
			if (ibits(fock_states(1,istate),isite-1,1) == 1) then         ! check if there is a up electron on that site
				do y=1,size(neighbours,2)                                 ! if so loop over all the sites that nieghbour it
					new_state(1) = IBCLR(fock_states(1,istate),isite-1)   ! remove the up electron from current site
					inbr = neighbours(isite,y)                            ! this the site number of the neighbour
					if (ibits(new_state(1),inbr-1,1) == 0) then           ! check if the neighbour site is empty (no up electron)
						new_state(2) = fock_states(2,istate)              ! the state after the hopping has same down component
						ne = 0                                            ! set the electron counter to zero (need to count them for anti-commutations)
						trans_site(1) = inbr; trans_site(2) = isite       ! the site it the elelectron started at and finished at
						do i=MINVAL(trans_site),MAXVAL(trans_site)-1
							ne = ne + ibits(new_state(1),i-1,1)           ! count the number of up electrons between start and finish site
							ne = ne + ibits(new_state(2),i-1,1)           ! count the number of down electrons between start and finish site
						end do
						if (MOD(ne,2) == 0) then                          ! if even amount of electrons between then anticommutation combine to 1
							phase = 1
						else
							phase = -1                                    ! if odd amount of electrons between then anticommutation combine to -1
						end if
						new_state(1) = ibset(new_state(1),inbr-1)         ! find the up component of the new fock state
						do j=1,nstates
							if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
								new_index = j                             ! search through the fock states to find the index of the new FS
							end if
						end do
						state_index = istate + 1 - mblock(n_up,n_dn)      ! find the row and column to add the 't' to the Hamiltonian matrix 
						new_index = new_index + 1 - mblock(n_up,n_dn)"

echo "						select case (n_up)"
for ((i=0; i<=nsites; i++)); do
	echo "							case($i)"
	echo "								select case (n_dn)"
	for ((j=0; j<=nsites; j++)); do
		echo "									case($j)"
		echo "										H$i$j(state_index,new_index) = t*phase"
	done
	echo "								end select"
done
echo -n "						end select"

echo "
					end if
				end do
			end if
			if (ibits(fock_states(2,istate),isite-1,1) == 1) then         ! Repeat identical process for hoppping of down electrons
				do y=1,size(neighbours,2)
					new_state(2) = IBCLR(fock_states(2,istate),isite-1)
					inbr = neighbours(isite,y)
					if (ibits(new_state(2),inbr-1,1) == 0) then
						new_state(1) = fock_states(1,istate)
						ne = 0
						trans_site(1) = inbr; trans_site(2) = isite
						do i=MINVAL(trans_site)+1,MAXVAL(trans_site)
							ne = ne + ibits(new_state(1),i-1,1)     ! count the number of up electrons in that state
							ne = ne + ibits(new_state(2),i-1,1)     ! count the number of down electrons in that state
						end do
						if (MOD(ne,2) == 0) then 
							phase = 1
						else
							phase = -1
						end if
						new_state(2) = ibset(new_state(2),inbr-1)
						do j=1,nstates
							if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
							new_index = j
							end if
						end do
						state_index = istate + 1 - mblock(n_up,n_dn)
						new_index = new_index + 1 - mblock(n_up,n_dn)"

echo "						select case (n_up)"
for ((i=0; i<=nsites; i++)); do
	echo "							case($i)"
	echo "								select case (n_dn)"
	for ((j=0; j<=nsites; j++)); do
		echo "									case($j)"
		echo "										H$i$j(state_index,new_index) = t*phase"
	done
	echo "								end select"
done
echo -n "						end select"

echo "
					end if
				end do
			end if
		end do
	end do
	end do
	end do

end subroutine make_hamiltonian2

!************************************************************************************

subroutine solve_hamiltonian1(E,U,mu)

	!  %--------------------------------------------------------------------------------%
	!  | This program makes the on diagonal terms and solves the lowest eigenvalue of   |
	!  | each submatrix.                                                                |
	!  |                                                                                |
	!  | The routines loops over all the submatrices and then over each state then      |
	!  | site. It counts the amount of up and down electrons on that site. It then adds |
	!  | the site potential of that site to the corresponding diagonal term in the      |
	!  | Hamiltonian matrix as many times as there are electrons. In addition, if there |
	!  | are two electrons (up and down) on the site an additional 'U' term is added.   |
	!  | The final step calls LAPACK to solve for the lowest eigenvalue of each         |
	!  | submatrix.                                                                     |
	!  %--------------------------------------------------------------------------------%

	implicit none

	integer :: n_up, n_dn,i,j,isite,istate                 ! counters for loops
	integer :: ne                                          ! counts the number of electrons
	real, intent(in) :: E(nsites)                          ! site potentials
	real, intent(in) :: mu                                 ! on-site coulomb interactions (defined at top of main.f90)
	real, intent(in) :: U                                  ! on-site coulomb interactions (defined at top of main.f90)
	real, allocatable, dimension(:) :: wtemp               ! temporary matrix to hold the eigenvalues
	real, allocatable, dimension(:,:) :: htemp, vtemp      ! temporary matrices to hold the hamiltonain and eigenvectors respectively

	do n_up=0,nsites
		do n_dn=n_up,nsites
			
			allocate(htemp(msize(n_up,n_dn),msize(n_up,n_dn)),wtemp(msize(n_up,n_dn)))  ! allocate temporary matrices of proper size
			allocate(vtemp(msize(n_up,n_dn),msize(n_up,n_dn)))
			"

echo "			select case (n_up)"
for ((i=0; i<=nsites; i++)); do
	echo "				case($i)"
	echo "					select case (n_dn)"
	for ((j=0; j<=nsites; j++)); do
		echo "						case($j)"
		echo "							htemp=H$i$j"
	done
	echo "					end select"
done
echo "			end select"

echo "
			wtemp=0.0
			do istate=1,msize(n_up,n_dn)              ! loop over all the states of the matrix
				do isite=1,nsites                     ! loop over all the site of each state
					ne=0
					ne = ne + IBITS(fock_states(1,istate+mblock(n_up,n_dn)-1),isite-1,1)  ! check if up electron on that site of that state
					ne = ne + IBITS(fock_states(2,istate+mblock(n_up,n_dn)-1),isite-1,1)  ! check if down electron on that site of that state
					htemp(istate,istate) = htemp(istate,istate) + ne*E(isite)             ! add ne*(the site potential of that site)
					if (ne == 2) then 
						htemp(istate,istate) = htemp(istate,istate) + U                   ! if there is both an up and down electron add a U
					end if
				end do
			end do
			call ssyevr_lapack1(msize(n_up,n_dn),htemp,wtemp,vtemp)                       ! call the LAPACK routine
			e_ground(n_up,n_dn) = wtemp(1) - mu*(n_up+n_dn)  ! grand potentials           ! record the smallest grand_potential (W(1) is smallest eigenvalue)
			deallocate(htemp,wtemp,vtemp)
		end do
	end do

end subroutine solve_hamiltonian1

!************************************************************************************

subroutine solve_hamiltonian2(E,U,mu,n_up,n_dn)

	!  %--------------------------------------------------------------------------------%
	!  |  This program is identical to solve_hamiltonian1 except it only solves one     |
	!  |  matrix instead of looping over all of them and calls a different lapack       |
	!  |  routine which solves all the eigenvalues and eigenvectors.                    |
	!  %--------------------------------------------------------------------------------%

	implicit none

	integer :: i,j,ne,isite,istate                   ! conters for loops
	real, intent(in) :: E(nsites)                    ! site potentials
	real, intent(in) :: mu                           ! chemical potential
	real, intent(in) :: U                            ! onsite coulomb interactions
	integer, intent(in) :: n_up,n_dn                 ! number of up and down electrons
	real, allocatable, dimension(:) :: wtemp
	real, allocatable, dimension(:,:) :: htemp, vtemp
			
	allocate(htemp(msize(n_up,n_dn),msize(n_up,n_dn)),wtemp(msize(n_up,n_dn)))
	allocate(vtemp(msize(n_up,n_dn),msize(n_up,n_dn)))
	"

echo "	select case (n_up)"
for ((i=0; i<=nsites; i++)); do
	echo "		case($i)"
	echo "			select case (n_dn)"
	for ((j=0; j<=nsites; j++)); do
		echo "				case($j)"
		echo "					htemp=H$i$j"
	done
	echo "			end select"
done
echo "	end select"

echo "
	wtemp=0.0
	do istate=1,msize(n_up,n_dn)
		do isite=1,nsites
			ne=0
			ne = ne + IBITS(fock_states(1,istate+mblock(n_up,n_dn)-1),isite-1,1)
			ne = ne + IBITS(fock_states(2,istate+mblock(n_up,n_dn)-1),isite-1,1)
			htemp(istate,istate) = htemp(istate,istate) + ne*E(isite)
			if (ne == 2) then 
				htemp(istate,istate) = htemp(istate,istate) + U
			end if
		end do
	end do
	call ssyevr_lapack(msize(n_up,n_dn),htemp,wtemp,vtemp)
	do i=1,msize(n_up,n_dn)
   		grand_potential(i+mblock(n_up,n_dn)-1) = wtemp(i) - mu*(n_up+n_dn)  ! grand potentials
	end do
	do i=1,msize(n_up,n_dn)
   		eigenvectors(i+mblock(n_up,n_dn)-1,1:msize(n_up,n_dn)) = vtemp(1:msize(n_up,n_dn),i)  ! eigenvectors
	end do
	deallocate(htemp,wtemp,vtemp)

end subroutine solve_hamiltonian2

!************************************************************************************

integer function choose(j,k)

    ! %---------------------------------------------%
	! |  simple choose function from statistics:    |
	! |                                             |
	! |                         i!                  |
	! |       choose(i,j) =  --------               |
	! |                      j!(i-j)!               |
	! %---------------------------------------------%

    implicit none
    
    integer :: j,k,i2
    integer (kind=8) :: tj,tk,tchoose

    tchoose = 1
    tj = j
    tk = k
    do i2 = tj-tk+1,tj
       tchoose = tchoose * i2
    end do
    do i2 = 2, tk
       tchoose = tchoose / i2
    end do
    choose = tchoose

end function choose

!************************************************************************************

end module"