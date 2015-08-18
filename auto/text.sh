#!/bin/bash

nsites=4
nstates="$(echo "$((4**$nsites))")"

echo "module routines

	use lapack

	implicit none

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

	!--------------Other Parameters and Variables-------------------------------
	integer, parameter :: nsites = $nsites             ! number of states
	integer, parameter :: nstates$nstates           ! size of the fock state (FS) basis
	real :: grand_potential_ground=0.0                 ! the lowest grand ensemble energy
	integer, dimension(2,nstates) :: fock_states       ! array that stores each FS represented in binary (see above)
	integer, dimension(0:2**nsites) :: states_order    ! orders the up part of FS and down part of FS before combining them
	integer :: ne                                      ! number of electrons (counts electrons in various routines)
	integer :: total_states_up                         ! total number of states with no down electrons
	integer, dimension(nsites,nstates) :: PES_down, PES_up   ! lookup tables for PES of many-body groundstate (MBG) transformations (ex. c|Psi0>)
	integer, dimension(nsites,nstates) :: IPES_down, IPES_up ! lookup tables for PES of MBG transformations (ex. c^{dagger}|Psi0>)
	integer, dimension(nsites,nstates) :: phase_PES_down, phase_PES_up    ! to get anticommutation sign right
	integer, dimension(nsites,nstates) :: phase_IPES_down, phase_IPES_up  ! to get anticommutation sign right

	!------------------Hamiltonian Submatrices---------------------------------
	"
max=0
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
	real, dimension(nstates,$max) :: eigenvectors          ! the many body eigenvectors (MBE) only coefficients of basis states with same n_up,n_dn as it
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

	character(len=50), intent(out) :: filename
	character(len=10) :: mu_str, t_str, U_str, W_str, s_str
	real, intent(in) :: t,U,mu,delta

	write(mu_str,'(F4.1)') mu
	write(t_str,'(I2)') nint(t)
	write(W_str,'(I2)') nint(delta) 
	write(U_str,'(I2)') nint(U)
	write(s_str,'(I1)') nsites

	write(filename,'(A)') trim(adjustl('data/a')) // trim(adjustl(s_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // '-dos+ipr_t' // trim(adjustl(t_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // 'U' // trim(adjustl(U_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // 'W' // trim(adjustl(W_str)) 
	write(filename,'(A)') trim(adjustl(filename)) // 'mu' // trim(adjustl(mu_str)) // '.dat'

end subroutine make_filename

!************************************************************************************

subroutine random_gen_seed()

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

  real, dimension(nsites), intent(out) :: E  ! site potentials
  real, intent(in) :: delta                  ! width of disorder
  real :: random                             ! random number
  integer :: i                               ! counter for loop

  do i=1,4
    call random_number(random)           ! gives a random number between 0 and 1
    E(i) = delta*(random - 0.5)          ! centers the random numbers about 0 and then expands width from 1 to delta
  end do

end subroutine site_potentials

!************************************************************************************

subroutine num_sites()

	!    %----------------------------------------------------------------------------------------------%
	!    | Creates and orders the fock state basis. As an example the order for 2 site system is:       |
	!    |                                                                                              |
	!    |                                                 i=                                           |
	!    |                   | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 	|
	!    |   ----------------|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|   |
	!    |   fock_states(1,i)| 0 | 0 | 0 | 0 | 1 | 2 | 1 | 2 | 1 | 2  | 1  | 2  | 3	 | 3  | 3  | 3  |	|
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

	total_states_up = 2**nsites     ! total number of fock states with n_dn = 0
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
 	! |	                     |0||   1   ||      2      |     3     | 4 |      |
 	! |                                                                       |
 	! %-----------------------------------------------------------------------%

	do istate=0,total_states_up-1
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
		do j=1,total_states_up
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
	!  %---------------------------------------------------------------------------------------------------%

	implicit none
	
	integer :: i,j,new_state(2), new_index, position, isite

	PES_up = 0; PES_down = 0
    IPES_up = 0; IPES_down = 0
    phase_PES_up = 0; phase_PES_down = 0
    phase_IPES_up = 0; phase_IPES_down = 0

	! make the PES_up matrices

	do position = 1,nsites
		do i=1,nstates
			ne = 0
			if (ibits(fock_states(1,i),position-1,1) == 1) then
				PES_up(position,i) = 0
				phase_PES_up(position,i) = 0 
			else
				do isite=position,nsites
					ne = ne + ibits(fock_states(1,i),isite-1,1)     ! count the number of up electrons in that state
					end do
				do isite=position,nsites
					ne = ne + ibits(fock_states(2,i),isite-1,1)     ! count the number of down electrons in that state
				end do
				if (MOD(ne,2) == 0) then 
					phase_PES_up(position,i) = 1
				else
					phase_PES_up(position,i) = -1
				end if
				new_state(1) = ibset(fock_states(1,i),position-1)
				new_state(2) = fock_states(2,i)
				do j=1,nstates
					if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
						new_index = j
					end if
				end do
				PES_up(position,i) = new_index
			end if
		end do
	end do

	do position = 1,nsites
		do i=1,nstates
			ne = 0
			if (ibits(fock_states(2,i),position-1,1) == 1) then
				PES_down(position,i) = 0
				phase_PES_down(position,i) = 0 
			else
				do isite=position+1,nsites
					ne = ne + ibits(fock_states(1,i),isite-1,1)     ! count the number of up electrons in that state
				end do
				do isite=position,nsites
					ne = ne + ibits(fock_states(2,i),isite-1,1)     ! count the number of down electrons in that state
				end do
				if (MOD(ne,2) == 0) then 
					phase_PES_down(position,i) = 1
				else
					phase_PES_down(position,i) = -1
				end if
				new_state(2) = ibset(fock_states(2,i),position-1)
				new_state(1) = fock_states(1,i)
				do j=1,nstates
					if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
						new_index = j
					end if
				end do
				PES_down(position,i) = new_index
			end if
		end do
	end do

	sites: do j=1,nsites  ! calculating the IPES matrices by making them the opposite of the PES
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

	! this program makes the hamiltonians off diagonal terms. On diagonal terms are added during each loop

	implicit none

	real, intent(in) :: t

	integer :: istate,isite, inbr,new_state(2),phase, ne, i, trans_site(2), new_index,j,state_index,y,n_up,n_dn
	"

for ((n_up=0;n_up<=nsites;n_up++)); do
echo -n "		"
for ((n_dn=0;n_dn<=nsites;n_dn++)); do
	echo -n "H$n_up$n_dn=0.0"
	if [[ n_dn != $nsites ]]; then
		echo -n "; "
	fi
done
echo ""
done

echo "
	call matrix_sizes()

	do n_up = 0,nsites
	do n_dn = 0,nsites
	do istate = mblock(n_up,n_dn),mblock(n_up,n_dn) + msize(n_up,n_dn)-1
		do isite = 1,nsites
			if (ibits(fock_states(1,istate),isite-1,1) == 1) then
				do y=1,size(neighbours,2)
					new_state(1) = IBCLR(fock_states(1,istate),isite-1)
					inbr = neighbours(isite,y)
					if (ibits(new_state(1),inbr-1,1) == 0) then
						new_state(2) = fock_states(2,istate)
						ne = 0
						trans_site(1) = inbr; trans_site(2) = isite
						do i=MINVAL(trans_site),MAXVAL(trans_site)-1
							ne = ne + ibits(new_state(1),i-1,1)     ! count the number of up electrons in that state
							ne = ne + ibits(new_state(2),i-1,1)     ! count the number of down electrons in that state
						end do
						if (MOD(ne,2) == 0) then 
							phase = 1
						else
							phase = -1
						end if
						new_state(1) = ibset(new_state(1),inbr-1)
						do j=1,nstates
							if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
							new_index = j
							end if
						end do
						state_index = istate + 1 - mblock(n_up,n_dn)
						new_index = new_index + 1 - mblock(n_up,n_dn)"

echo "							select case (n_up)"
for ((i=0; i<=nsites; i++)); do
echo "								case($i)"
echo "									select case (n_dn)"
for ((j=0; j<=nsites; j++)); do
	echo "										case($j)"
	echo "											H$i$j(state_index,new_index) = t*phase"
done
echo "									end select"
done
echo -n "							end select"

echo "
					end if
				end do
			end if
			if (ibits(fock_states(2,istate),isite-1,1) == 1) then
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

echo "							select case (n_up)"
for ((i=0; i<=nsites; i++)); do
echo "								case($i)"
echo "									select case (n_dn)"
for ((j=0; j<=nsites; j++)); do
	echo "										case($j)"
	echo "											H$i$j(state_index,new_index) = t*phase"
done
echo "									end select"
done
echo -n "							end select"

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

	! this program makes the on diagonal terms and solves teh eigenenergies

	implicit none

	integer :: n_up, n_dn,i,j,ne,isite,istate

	real, intent(in) :: E(nsites)
		real, intent(in) :: mu 
		real, intent(in) :: U
		real, allocatable, dimension(:) :: wtemp
		real, allocatable, dimension(:,:) :: htemp, vtemp

	do n_up=0,nsites
		do n_dn=0,nsites
			
			allocate(htemp(msize(n_up,n_dn),msize(n_up,n_dn)),wtemp(msize(n_up,n_dn)))
			allocate(vtemp(msize(n_up,n_dn),msize(n_up,n_dn)))
			"

echo "				select case (n_up)"
for ((i=0; i<=nsites; i++)); do
echo "					case($i)"
echo "						select case (n_dn)"
for ((j=0; j<=nsites; j++)); do
	echo "							case($j)"
	echo "								htemp=H$i$j"
done
echo "						end select"
done
echo "				end select"

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

			call ssyevr_lapack1(msize(n_up,n_dn),htemp,wtemp,vtemp)

			e_ground(n_up,n_dn) = wtemp(1) - mu*(n_up+n_dn)  ! grand potentials

			deallocate(htemp,wtemp,vtemp)

		end do
	end do

end subroutine solve_hamiltonian1

!************************************************************************************

subroutine solve_hamiltonian2(E,U,mu,n_up,n_dn)

	! this program makes the on diagonal terms and solves it

	implicit none

	integer :: i,j,ne,isite,istate

	real, intent(in) :: E(nsites)
		real, intent(in) :: mu 
		real, intent(in) :: U
		integer, intent(in) :: n_up,n_dn
		real, allocatable, dimension(:) :: wtemp
		real, allocatable, dimension(:,:) :: htemp, vtemp
			
	allocate(htemp(msize(n_up,n_dn),msize(n_up,n_dn)),wtemp(msize(n_up,n_dn)))
	allocate(vtemp(msize(n_up,n_dn),msize(n_up,n_dn)))
	"

echo "		select case (n_up)"
for ((i=0; i<=nsites; i++)); do
echo "			case($i)"
echo "				select case (n_dn)"
for ((j=0; j<=nsites; j++)); do
	echo "					case($j)"
	echo "						htemp=H$i$j"
done
echo "				end select"
done
echo "		end select"

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
