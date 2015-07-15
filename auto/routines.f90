module routines

	implicit none

	integer, parameter :: nsites = 2
	integer, parameter :: int_kind = 4
	integer, parameter :: real_kind = kind(1.d0)
	integer, parameter :: STDOUT = 10
	real (real_kind), parameter :: PI_ = 4d0*atan(1d0)

	integer, allocatable, dimension(:,:) :: fock_states
	integer, allocatable, dimension(:) :: states_order
	integer :: ne
	integer :: total_states_up
	integer, dimension(nsites,4**nsites) :: PES_down=0, PES_up=0, IPES_down=0, IPES_up=0  !matrices for PES and IPES 
	integer, dimension(nsites,4**nsites) :: phase_PES_down=0, phase_PES_up=0, phase_IPES_down=0, phase_IPES_up=0  !do get anticommutation sign right

	integer, allocatable, dimension(:) :: block, temp_block, nstates_up

contains

	!------------------------------------------------------

	subroutine num_sites()

		implicit none

		integer :: i,j, istate, isite

		integer :: max_electrons

		total_states_up = 2**nsites
		max_electrons = nsites

		allocate(states_order(total_states_up))
		allocate(block(0:max_electrons),temp_block(0:max_electrons))
		allocate(nstates_up(0:nsites))

		block = 0  
	    nstates_up = 0
	    block(0) = 1
	    nstates_up(0) = 1

		do ne = 1, max_electrons
       		nstates_up(ne) = choose(nsites,ne)          ! number of spin-up fock states with ne electrons
       		block(ne) = block(ne-1) + nstates_up(ne-1) ! block index updating
    	end do

    	temp_block = block

		do istate=0,total_states_up-1
			ne = 0
			do isite=1,nsites
				ne = ne + ibits(istate,isite-1,1)     ! count the number of electrons in that state
			end do 
			states_order(temp_block(ne)) = istate
			temp_block(ne) = temp_block(ne) + 1
		end do

		deallocate(temp_block)

		allocate(fock_states(2,4**nsites))

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

	!------------------------------------------------------

	subroutine transformations()

		implicit none
		
		integer :: i,j,new_state(2), new_index, position, isite

		! make the PES_up matrices

		do position = 1,nsites
			do i=1,4**nsites
				ne = 0
				if (ibits(fock_states(1,i),position-1,1) == 1) then
					PES_up(position,i) = 0
					phase_PES_up(position,i) = 0 
				else
					do isite=position,nsites
						ne = ne + ibits(fock_states(1,isite),isite-1,1)     ! count the number of up electrons in that state
						ne = ne + ibits(fock_states(2,isite),isite-1,1)     ! count the number of down electrons in that state
					end do
					if (MOD(ne,2) == 1) then 
						phase_PES_up(position,i) = -1
					else
						phase_PES_up(position,i) = 1
					end if
					new_state(1) = ibset(fock_states(1,i),position-1)
					new_state(2) = fock_states(2,i)
					do j=1,4**nsites
						if(fock_states(1,j) == new_state(1) .and. fock_states(2,j) == new_state(2)) then
							new_index = j
						end if
					end do
					PES_up(position,i) = new_index
				end if
			end do
		end do

	end subroutine transformations

	!------------------------------------------------------

	subroutine new_state(old,new)

		implicit none
		
		integer, intent(in) :: old
		integer, intent(out) :: new 




	end subroutine new_state

	!-------------------------------------------------------

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

	!------------------------------------------------------

end module